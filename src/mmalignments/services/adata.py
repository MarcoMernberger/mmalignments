import numpy as np
import pandas as pd
import anndata as ad
from anndata import AnnData
from scipy.sparse import csr_matrix
from pandas import DataFrame
from mmalignments.services.io import from_json, parents
from pathlib import Path
from mmalignments.models.tags import ElementTag, Method, PartialElementTag, Stage, State
from mmalignments.models.elements import Element, AdataElement, generate_element_key_name, element
from typing import Iterable, Callable
from mmalignments import __version__


def get_columns_for_roles(column_roles: dict[str, str], role: str) -> list[str]:
    return [col for col, r in column_roles.items() if r == role]


def create_adata_from_frame(
    outfile: str | Path,
    data_table: str | Path,
    samples_table: str | Path,
    roles_json: str | Path,
    *,
    genes_table: str | Path | None,
    layers: tuple[str, ...] = ("raw", "cpm", "vst"),
    default_layer: str = "vst",
    gene_metadata_role: str | None = "annotation",
    make_sparse: bool = False,
    samples_index: str = "sample",
    genes_index: str = "gene_stable_id",
    dtype: str = "float32",
    rename_columns: Callable | dict[str, Callable] | None = None,
):
    def __write():
        print("run is actually called")
        anndata = build_from_files(
            genes_table=genes_table,
            samples_table=samples_table,
            data_table=data_table,
            roles_json=roles_json,
            layers=layers,
            default_layer=default_layer,
            gene_metadata_role=gene_metadata_role,
            make_sparse=make_sparse,
            samples_index=samples_index,
            genes_index=genes_index,
            dtype=dtype,
            rename_columns=rename_columns,
        )
        print("anndata", anndata)
        parents(outfile)
        anndata.write_h5ad(filename=outfile, convert_strings_to_categoricals=True)
    
    return __write


def build_from_files(
        genes_table: str | Path | None,
        samples_table: str | Path,
        data_table: str | Path,
        roles_json: str | Path,
        layers: tuple[str, ...] = ("raw", "cpm", "vst"),
        default_layer: str = "vst",
        gene_metadata_role: str | None = "annotation",
        make_sparse: bool = False,
        samples_index: str = "sample",
        genes_index: str = "gene_stable_id",
        dtype: str = "float32",
        rename_columns: Callable | dict[str, Callable] | None = None,
):
    def default_rename_fn(c: str) -> str:
        return c.split()[0]
    
    def __build():
        rename_fn = rename_columns or default_rename_fn
        genes = pd.read_csv(genes_table, sep="\t") if genes_table else None
        samples = pd.read_csv(samples_table, sep="\t")
        samples.index = samples[samples_index].astype(str)  # separate index
        data_in = pd.read_csv(data_table, sep="\t")
        data_in.index = data_in[genes_index].astype(str)    # separate index
        roles = from_json(Path(roles_json))
        data = {}
        print(samples)
        for role in layers:
            columns_for_role = get_columns_for_roles(roles, role)
            data_view = data_in[columns_for_role]
            data_view = data_view.rename(columns=rename_fn)
            data[role] = data_view
        if not genes:
            gene_metadata_columns = get_columns_for_roles(roles, gene_metadata_role) if gene_metadata_role else None

        return build_anndata(
            data=data,
            data_obs=samples,
            data_var=genes,
            layers=layers,
            default_layer=default_layer,
            gene_metadata_columns=gene_metadata_columns,
            make_sparse=make_sparse,
            dtype=dtype,
        )
    return __build()


def build_anndata(
    data,
    data_obs=None,
    data_var=None,
    *,
    layers=("raw", "cpm", "vst"),
    default_layer="raw",
    gene_metadata_columns=None,
    make_sparse=False,
    dtype="float32",
) -> AnnData:
    """
    Builds an AnnData object from DataFrames.

    Expected structure
    ------------------
    data:
        Either:
        1. dict[str, pd.DataFrame]
           Example:
               {
                   "raw": raw_df,
                   "cpm": cpm_df,
                   "vst": vst_df
               }

           Each matrix:
               rows    = Genes
               columns = Samples

        OR:

        2. A single pd.DataFrame with MultiIndex columns:
               level 0 = layer/view
               level 1 = sample

    data_obs:
        Sample metadata.
        Index must contain sample IDs.

    data_var:
        Gene metadata.
        Index must contain gene IDs.

    gene_metadata_columns:
        If data_var=None:
        List of column names from `data`,
        to be extracted as gene metadata.

    make_sparse:
        Whether matrices should be stored as CSR sparse.

    dtype:
        Data type of the matrix.
    """

    # ------------------------------------------------------------
    # Normalize input
    # ------------------------------------------------------------

    if isinstance(data, dict):

        matrices = data

    elif isinstance(data, pd.DataFrame):

        if not isinstance(data.columns, pd.MultiIndex):
            raise ValueError(
                "If `data` is a DataFrame, "
                "the columns must be a MultiIndex: "
                "(layer, sample)"
            )

        matrices = {
            layer: data[layer]
            for layer in data.columns.get_level_values(0).unique()
        }

    else:
        raise TypeError(
            "`data` must be a dict or DataFrame."
        )

    # ------------------------------------------------------------
    # Validate layers
    # ------------------------------------------------------------

    available_layers = list(matrices.keys())

    for layer in layers:
        if layer not in matrices:
            raise ValueError(
                f"Layer '{layer}' not found. "
                f"Available: {available_layers}"
            )

    if default_layer not in matrices:
        raise ValueError(
            f"default_layer '{default_layer}' not found."
        )

    # ------------------------------------------------------------
    # Reference matrix
    # ------------------------------------------------------------
    ref = matrices[default_layer]

    # ------------------------------------------------------------
    # Extract gene metadata from data
    # ------------------------------------------------------------

    if data_var is None and gene_metadata_columns is not None:

        existing_cols = [
            c for c in gene_metadata_columns
            if c in ref.columns
        ]

        if existing_cols:

            data_var = (
                ref[existing_cols]
                .copy()
            )

            ref = ref.drop(columns=existing_cols)

            for layer in matrices:
                matrices[layer] = matrices[layer].drop(
                    columns=existing_cols,
                    errors="ignore"
                )

    # ------------------------------------------------------------
    # Shapes prüfen
    # ------------------------------------------------------------

    genes = ref.index
    samples = ref.columns

    for layer_name, df in matrices.items():

        if not genes.equals(df.index):
            raise ValueError(
                f"Genes do not match in layer '{layer_name}'."
            )

        if not samples.equals(df.columns):
            print("Expected samples:", samples)
            print("Found samples:", df.columns)
            raise ValueError(
                f"Samples do not match in layer '{layer_name}'."
            )

    # ------------------------------------------------------------
    # Create obs
    # ------------------------------------------------------------

    if data_obs is None:

        data_obs = pd.DataFrame(index=samples)

    else:

        data_obs = data_obs.copy()

        missing = samples.difference(data_obs.index)

        if len(missing) > 0:
            raise ValueError(
                f"Samples missing in data_obs: {list(missing)}"
            )

        data_obs = data_obs.loc[samples]

    # ------------------------------------------------------------
    # Create var
    # ------------------------------------------------------------

    if data_var is None:

        data_var = pd.DataFrame(index=genes)

    else:

        data_var = data_var.copy()

        missing = genes.difference(data_var.index)

        if len(missing) > 0:
            raise ValueError(
                f"Genes missing in data_var: {list(missing)}"
            )

        data_var = data_var.loc[genes]

    # ------------------------------------------------------------
    # Optimize categories
    # ------------------------------------------------------------

    for col in data_obs.columns:

        if data_obs[col].dtype == "object":

            nunique = data_obs[col].nunique()

            if nunique < len(data_obs) * 0.5:

                data_obs[col] = (
                    data_obs[col]
                    .astype("category")
                )

    # ------------------------------------------------------------
    # Prepare main matrix
    # ------------------------------------------------------------

    X = ref.T.values.astype(dtype)

    if make_sparse:
        X = csr_matrix(X)

    # ------------------------------------------------------------
    # Create AnnData
    # ------------------------------------------------------------

    adata = ad.AnnData(
        X=X,
        obs=data_obs,
        var=data_var,
    )

    # ------------------------------------------------------------
    # Add layers
    # ------------------------------------------------------------

    for layer_name, df in matrices.items():

        arr = df.T.values.astype(dtype)

        if make_sparse:
            arr = csr_matrix(arr)

        adata.layers[layer_name] = arr

    # ------------------------------------------------------------
    # Document default layer
    # ------------------------------------------------------------

    adata.uns["default_layer"] = default_layer

    return adata


@element
def AdataFromFrame(
    data_tsv: str | Path,
    samples_tsv: str | Path,
    roles_json: str | Path,
    *,
    tag: PartialElementTag | ElementTag | None = None,
    outdir: str | Path | None = None,
    filename: str | None = None,
    root: str | None = None,
    genes_tsv: str | Path | None = None,
    layers: tuple[str, ...] = ("raw", "cpm"),
    default_layer: str = "raw",
    gene_metadata_role: str | None = "annotation",
    # gene_roles: str | dict[str, str] | None  = "annotation",
    sample_roles: str | dict[str, str] | None = None,
    make_sparse: bool = False,
    dtype: str = "float32",
):
    root = root or Path(data_tsv).stem

    tag = ElementTag(
        root=root,
        level=0,
        omics=None,
        stage=Stage.INPUT,
        method=Method.CHECK,
        state=State.RAW,
        ext=".tsv",
    ).merge(tag)
    key, name = generate_element_key_name(
        tag=tag,
        tool_name="mmalignments",
        tool_version=__version__,
    )
    inputs=(data_tsv, samples_tsv, roles_json)
    if genes_tsv:
        inputs += (genes_tsv,)

    determinants = layers + (default_layer, str(make_sparse), str(dtype))  #  + (gene_roles, ) if isinstance(gene_roles, str) else tuple(gene_roles or ())
    outdir = outdir or Path("cache") / "adata"
    filename = filename or f"{key}.h5ad"
    h5ad = outdir / filename
    runner = create_adata_from_frame(
        h5ad,
        genes_table=genes_tsv,
        samples_table=samples_tsv,
        data_table=data_tsv,
        roles_json= roles_json,
        layers=layers,
        default_layer=default_layer,
        gene_metadata_role=gene_metadata_role,
    )
    return AdataElement(
        key=key,
        run=runner,
        tag=tag,
        root=root,
        h5ad=h5ad,
        # obs_roles=sample_roles,
        # var_roles=gene_roles,
        determinants=determinants,
        inputs=inputs,
        pres=None,
        name=name,
    )
