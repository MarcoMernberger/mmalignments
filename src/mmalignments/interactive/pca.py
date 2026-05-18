import numpy as np
import pandas as pd
import panel as pn
import param

from datetime import datetime
from pathlib import Path

from sklearn.decomposition import PCA

import plotly.express as px


class PCAAnalysis(param.Parameterized):


    # ==============================================================
    # PCA parameters
    # ==============================================================

    pca_name = param.String(default="default_pca")

    layer = param.Selector(
        default=None,
        objects=[],
    )

    n_components = param.Integer(
        default=5,
        bounds=(2, 50),
    )

    pc_x = param.Integer(
        default=1,
        bounds=(1, 50),
    )

    pc_y = param.Integer(
        default=2,
        bounds=(1, 50),
    )

    # ==============================================================
    # Gene selection
    # ==============================================================

    gene_subset_mode = param.Selector(
        default="all",
        objects=[
            "all",
            "manual",
            "top_variable",
        ],
    )

    top_variable_genes = param.Integer(
        default=500,
        bounds=(10, 50000),
    )

    gene_selection = param.ListSelector(
        default=[],
        objects=[],
    )

    # ==============================================================
    # Preprocessing
    # ==============================================================

    center = param.Boolean(default=True)

    scale = param.Boolean(default=False)

    # ==============================================================
    # Coloring
    # ==============================================================

    color_by = param.Selector(
        default=None,
        objects=[],
    )

    # ==============================================================
    # Labels
    # ==============================================================

    show_labels = param.Boolean(default=False)

    label_size = param.Integer(
        default=10,
        bounds=(6, 40),
    )

    label_distance = param.Number(
        default=0.0,
        bounds=(0.0, 5.0),
    )

    # ==============================================================
    # Runtime
    # ==============================================================

    pca_updated = param.Event()

    auto_recompute = param.Boolean(default=True)

    # ==============================================================
    # Saving
    # ==============================================================

    output_dir = param.String(default="results")

    output_filename = param.String(default="analysis.h5ad")

    save_status = param.String(default="")

    save_fig_status = param.String(default="")

    # ==============================================================
    # Init
    # ==============================================================

    def __init__(self, adata, **params):

        super().__init__(**params)

        self.adata = adata

        self.user_defined_pca_name = (
            self.pca_name != "default_pca"
        )

        # ----------------------------------------------------------
        # layers
        # ----------------------------------------------------------

        layers = list(adata.layers.keys())

        if adata.X is not None:
            layers = ["X"] + layers

        self.param.layer.objects = layers

        if self.layer is None and layers:
            self.layer = layers[0]

        # ----------------------------------------------------------
        # obs columns
        # ----------------------------------------------------------

        obs_columns = list(adata.obs.columns)

        self.param.color_by.objects = obs_columns

        if obs_columns and self.color_by is None:
            self.color_by = obs_columns[0]

        # ----------------------------------------------------------
        # genes
        # ----------------------------------------------------------

        self.param.gene_selection.objects = (
            list(adata.var_names)
        )

        # ----------------------------------------------------------
        # watches
        # ----------------------------------------------------------

        self.param.watch(
            self._parameter_changed,
            [
                "layer",
                "n_components",
                "gene_subset_mode",
                "top_variable_genes",
                "gene_selection",
                "center",
                "scale",
            ],
        )

        # ----------------------------------------------------------
        # initial PCA
        # ----------------------------------------------------------

        self.compute_pca()

    # ==============================================================
    # Dynamic bounds
    # ==============================================================

    @param.depends("n_components", watch=True)
    def _update_component_bounds(self):

        bounds = (1, self.n_components)

        self.param.pc_x.bounds = bounds
        self.param.pc_y.bounds = bounds

        self.pc_x = min(
            self.pc_x,
            self.n_components,
        )

        self.pc_y = min(
            self.pc_y,
            self.n_components,
        )

    # ==============================================================
    # Name tracking
    # ==============================================================

    @param.depends("pca_name", watch=True)
    def _track_custom_name(self):

        if self.pca_name != self.automatic_pca_name:

            self.user_defined_pca_name = True

    # ==============================================================
    # Helpers
    # ==============================================================

    @property
    def output_path(self):

        return (
            Path(self.output_dir)
            / self.output_filename
        )

    @property
    def automatic_pca_name(self):

        parts = [
            str(self.layer),
            str(self.gene_subset_mode),
        ]

        if self.gene_subset_mode == "top_variable":

            parts.append(
                f"hvg{self.top_variable_genes}"
            )

        if self.scale:
            parts.append("scaled")

        if self.center:
            parts.append("centered")

        return "_".join(parts)

    @property
    def current_embedding_key(self):

        return f"X_pca_{self.pca_name}"

    @property
    def current_loading_key(self):

        return f"PCs_{self.pca_name}"

    @property
    def current_uns_key(self):

        return f"pca_{self.pca_name}"

    # ==============================================================
    # Parameter change callback
    # ==============================================================

    def _parameter_changed(self, event):

        if self.auto_recompute:

            if not self.user_defined_pca_name:

                self.pca_name = (
                    self.automatic_pca_name
                )

            self.compute_pca()

    # ==============================================================
    # Matrix extraction
    # ==============================================================

    def get_matrix(self):

        if self.layer == "X":

            X = self.adata.X

        else:

            X = self.adata.layers[self.layer]

        if hasattr(X, "toarray"):

            X = X.toarray()

        genes = self.select_genes(X)

        return X[:, genes], genes

    # ==============================================================
    # Gene selection
    # ==============================================================

    def select_genes(self, X):

        if self.gene_subset_mode == "all":

            return np.arange(X.shape[1])

        if self.gene_subset_mode == "manual":

            if not self.gene_selection:

                return np.arange(X.shape[1])

            mask = self.adata.var_names.isin(
                self.gene_selection
            )

            return np.where(mask)[0]

        if self.gene_subset_mode == "top_variable":

            variances = np.var(X, axis=0)

            top_idx = np.argsort(
                variances
            )[-self.top_variable_genes:]

            return np.sort(top_idx)

        return np.arange(X.shape[1])

    # ==============================================================
    # PCA computation
    # ==============================================================

    def compute_pca(self):

        X, gene_idx = self.get_matrix()

        # ----------------------------------------------------------
        # center
        # ----------------------------------------------------------

        if self.center:

            X = X - X.mean(axis=0)

        # ----------------------------------------------------------
        # scale
        # ----------------------------------------------------------

        if self.scale:

            std = X.std(axis=0)

            std[std == 0] = 1

            X = X / std

        # ----------------------------------------------------------
        # safe n_components
        # ----------------------------------------------------------

        max_components = min(
            X.shape[0],
            X.shape[1],
        )

        n_components = min(
            self.n_components,
            max_components,
        )

        self.n_components = n_components

        # ----------------------------------------------------------
        # PCA
        # ----------------------------------------------------------

        pca = PCA(
            n_components=n_components
        )

        coords = pca.fit_transform(X)

        # ----------------------------------------------------------
        # automatic name
        # ----------------------------------------------------------

        if not self.user_defined_pca_name:

            self.pca_name = (
                self.automatic_pca_name
            )

        # ----------------------------------------------------------
        # keys
        # ----------------------------------------------------------

        embedding_key = (
            self.current_embedding_key
        )

        loading_key = (
            self.current_loading_key
        )

        uns_key = (
            self.current_uns_key
        )

        # ----------------------------------------------------------
        # store coordinates
        # ----------------------------------------------------------

        self.adata.obsm[
            embedding_key
        ] = coords

        # ----------------------------------------------------------
        # loadings
        # ----------------------------------------------------------

        loadings = np.full(
            (
                self.adata.n_vars,
                n_components,
            ),
            np.nan,
        )

        loadings[
            gene_idx,
            :
        ] = pca.components_.T

        self.adata.varm[
            loading_key
        ] = loadings

        # ----------------------------------------------------------
        # metadata
        # ----------------------------------------------------------

        metadata = {

            "variance_ratio": (
                pca.explained_variance_ratio_
            ),

            "variance": (
                pca.explained_variance_
            ),

            "params": {

                "layer": self.layer,

                "n_components": (
                    n_components
                ),

                "gene_subset_mode": (
                    self.gene_subset_mode
                ),

                "top_variable_genes": (
                    self.top_variable_genes
                ),

                "center": self.center,

                "scale": self.scale,
            },

            "gene_indices": (
                gene_idx.tolist()
            ),

            "created": (
                datetime.now().isoformat()
            ),
        }

        self.adata.uns[
            uns_key
        ] = metadata

        # ----------------------------------------------------------
        # registry
        # ----------------------------------------------------------

        if "pca_registry" not in self.adata.uns:

            self.adata.uns[
                "pca_registry"
            ] = {}

        self.adata.uns[
            "pca_registry"
        ][self.pca_name] = {

            "obsm": embedding_key,

            "varm": loading_key,

            "uns": uns_key,

            "created": (
                datetime.now().isoformat()
            ),

            "params": metadata["params"],
        }

        # ----------------------------------------------------------
        # active PCA
        # ----------------------------------------------------------

        self.adata.uns[
            "active_pca"
        ] = self.pca_name

        # ----------------------------------------------------------
        # trigger UI
        # ----------------------------------------------------------

        self.param.trigger(
            "pca_updated"
        )

    # ==============================================================
    # PCA dataframe
    # ==============================================================

    def pca_dataframe(self):

        coords = self.adata.obsm[
            self.current_embedding_key
        ]

        explained = self.adata.uns[
            self.current_uns_key
        ]["variance_ratio"]
        x = coords[:, self.pc_x - 1]
        y = coords[:, self.pc_y - 1]
        angles = np.linspace(
            0,
            2 * np.pi,
            len(x),
            endpoint=False,
        )

        label_x = x + (
            np.cos(angles)
            * self.label_distance
        )

        label_y = y + (
            np.sin(angles)
            * self.label_distance
        )

        df = pd.DataFrame({

            "PCX": x,

            "PCY": y,

            "sample": (
                self.adata.obs_names
            ),

            "label_x": label_x,

            "label_y": label_y,
        })

        if self.color_by is not None:

            df["color"] = (
                self.adata.obs[
                    self.color_by
                ]
                .astype(str)
                .values
            )

        df.attrs["explained"] = explained

        return df

    # ==============================================================
    # Plot
    # ==============================================================

    def _build_fig(self, df):

        explained = df.attrs["explained"]

        fig = px.scatter(
            df,
            x="PCX",
            y="PCY",
            color="color" if "color" in df else None,
            hover_name="sample",
            title=f"PCA ({self.layer})",
        )

        if self.show_labels:
            fig.add_scatter(
                x=df["label_x"],
                y=df["label_y"],
                mode="text",
                text=df["sample"],
                showlegend=False,
            )

        fig.update_layout(height=700, template="plotly_white")

        return fig

    @param.depends(
        "pc_x",
        "pc_y",
        "color_by",
        "show_labels",
        "label_size",
        "label_distance",
        "pca_updated",
    )
    def plot(self):

        df = self.pca_dataframe()

        explained = df.attrs[
            "explained"
        ]

        fig = px.scatter(
            df,
            x="PCX",
            y="PCY",
            color=(
                "color"
                if "color" in df.columns
                else None
            ),
            hover_name="sample",
            title=f"PCA ({self.pca_name})",
            labels={

                "PCX": (
                    f"PC{self.pc_x} "
                    f"({explained[self.pc_x - 1]:.1%})"
                ),

                "PCY": (
                    f"PC{self.pc_y} "
                    f"({explained[self.pc_y - 1]:.1%})"
                ),
            },
        )

        if self.show_labels:

            fig.add_scatter(
                x=df["label_x"],
                y=df["label_y"],
                mode="text",
                text=df["sample"],
                showlegend=False,
            )
            for _, row in df.iterrows():

                fig.add_annotation(

                    x=row["PCX"],
                    y=row["PCY"],

                    ax=row["label_x"],
                    ay=row["label_y"],

                    xref="x",
                    yref="y",

                    axref="x",
                    ayref="y",

                    showarrow=True,

                    arrowhead=1,

                    arrowsize=1,

                    arrowwidth=1,

                    text="",
                )

        fig.update_layout(
            height=700,
            template="plotly_white",
        )

        fig.update_traces(
            textposition="top center",
            textfont_size=self.label_size,
        )

        return fig

    @param.depends("pca_updated")
    def scree_plot(self):

        key = f"pca_{self.pca_name}"
        info = self.adata.uns.get(key, None)

        if info is None:
            return None

        evr = np.array(info["explained_variance_ratio"])
        cum = np.cumsum(evr)

        x = np.arange(1, len(evr) + 1)

        fig = px.line(
            x=x,
            y=evr,
            markers=True,
            labels={"x": "PC", "y": "Explained variance"},
            title="Scree Plot (Explained Variance)"
        )

        fig.add_scatter(
            x=x,
            y=cum,
            mode="lines",
            name="Cumulative"
        )

        fig.update_layout(height=350)

        return fig
    # ==============================================================
    # Saving
    # ==============================================================

    def save(self, path):

        self.adata.write_h5ad(path)

    def save_adata(self, event=None):

        self.output_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.adata.write_h5ad(
            self.output_path
        )

        self.save_status = (
            f"Saved to {self.output_path}"
        )

    def save_figure(self, event=None):

        df = self.pca_dataframe()

        fig = self._build_fig(df)   # wir extrahieren Plot-Logik

        base = self.output_path.with_suffix("")
        base = Path(base)

        fig.write_image(f"{base}.png")
        fig.write_image(f"{base}.svg")
        fig.write_image(f"{base}.pdf")

        self.save_fig_status = f"Figure saved to {base}"
    # ==============================================================
    # Panel
    # ==============================================================

    def panel(self):

        save_button = pn.widgets.Button(
            name="Save AnnData",
            button_type="primary",
        )
        save_button.on_click(
            self.save_adata
        )
        save_fig_button = pn.widgets.Button(
            name="Save Figure",
            button_type="success",
        )
        save_fig_button.on_click(self.save_figure)

        controls = pn.Param(
            self,
            widgets={

                "gene_selection": {

                    "type": (
                        pn.widgets.MultiChoice
                    ),

                    "solid": False,
                }
            },
            sizing_mode="stretch_width",
        )

        save_status = pn.pane.Markdown(self.param.save_status)
        save_fig_status = pn.pane.Markdown(self.param.save_fig_status)

        # plot = pn.panel(self.plot)

        return pn.template.FastListTemplate(

            title="PCA Explorer",

            sidebar=[

                controls,

                pn.pane.Markdown(
                    "## Save"
                ),

                save_button,

                save_status,

                pn.pane.Markdown("## Export"),
                
                save_fig_button,
                
                save_fig_status,
            ],

            main=[
                pn.Row(
                    pn.panel(self.plot),
                    pn.panel(self.scree_plot),
                )
            ],
        )