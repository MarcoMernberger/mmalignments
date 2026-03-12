from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterator, Mapping

from mmalignments.services.io import load_param_json
from mmalignments.models.resources import ResourceConfig  # type: ignore[import]


###############################################################################
# Params
###############################################################################


@dataclass(frozen=True)
class Params:
    _override_params: Mapping[str, Any] = field(default_factory=dict)

    def __init__(self, **kwargs: Any) -> None:
        object.__setattr__(self, "_override_params", dict(kwargs))

    @classmethod
    def of(cls, **kwargs: Any) -> "Params":
        return cls(**kwargs)

    def __getitem__(self, key: str) -> Any:
        return self._override_params[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._override_params)

    def __len__(self) -> int:
        return len(self._override_params)

    def to_dict(self, *, drop_none: bool = True) -> dict[str, Any]:
        if not drop_none:
            return dict(self._override_params)
        return {k: v for k, v in self._override_params.items() if v is not None}

    def get(self, key: str, default: Any = None) -> Any:
        return self._override_params.get(key, default)

    def items(self) -> Iterator[tuple[str, Any]]:
        return iter(self._override_params.items())

    def __contains__(self, key: str) -> bool:
        return key in self._override_params

    def override(self, **kwargs: Any) -> "Params":
        return Params(**{**self._override_params, **kwargs})

    def __repr__(self) -> str:
        return f"Params({self._override_params})"

    def __str__(self) -> str:
        return self.__repr__()


###############################################################################
# Rendering
###############################################################################
RenderFn = Callable[[str | None, Any], list[str]]


def render_flag(flag: str | None, value: Any) -> list[str]:
    return [flag] if value is True and flag else []


def render_value(flag: str | None, value: Any) -> list[str]:
    if flag is None:
        return [str(value)]
    return [flag, str(value)]


def render_multi_repeat(flag: str | None, value: list[Any]) -> list[str]:
    if flag is None:
        return [str(v) for v in value]
    out: list[str] = []
    for v in value:
        out.extend([flag, str(v)])
    return out


_RENDER: dict[str, Callable] = {
    "value": render_value,
    "flag": render_flag,
    "multi_repeat": render_multi_repeat,
}


def _parse_render(r: Any) -> Callable:
    if r is None:
        return render_value
    if isinstance(r, str):
        if r not in _RENDER:
            raise ValueError(f"Unknown render '{r}'. Allowed: {sorted(_RENDER)}")
        return _RENDER[r]
    raise ValueError(f"Invalid render format: {r!r}")


###############################################################################
# ParamSpec and ParamSet
###############################################################################


@dataclass(frozen=True)
class ParamSpec:
    name: str
    flag: str | None  # "--threads" or None for positional
    dtype: type | tuple[type, ...]  # expected type(s)
    required: bool = False  # is this required?
    default: Any = None  # default value
    affects_output: bool = True  # do we need a re run?
    description: str = ""  # for documentation
    render: Callable = render_value  # "flag", "value", "multi"


#    normalize: Callable[[Any], Any] | None = None  # optional normalization function
###############################################################################
# ToolThreadSpec
###############################################################################


@dataclass(frozen=True)
class ToolThreadSpec(ParamSpec):
    """Declares how a specific tool or subroutine uses threads.

    This is the **tool-side** half of thread allocation. It lives on the
    External subclass and describes *capability*, not a concrete value.

    Parameters
    ----------
    param_flag : str | None
        CLI flag to pass the thread count to the tool (e.g. ``"-t"``,
        ``"--threads"``, ``"--nthreads"``). ``None`` means the tool does not
        accept a thread flag but still benefits from CPU reservation.
    max_threads : int | None
        Hard cap on threads for this tool/subroutine. ``None`` = no cap,
        use all available.
    fraction : float
        Fraction of available CPUs to use (default 1.0 = all available).
        E.g. ``0.5`` means "use half the machine's CPUs".
    multi : bool
        Whether the tool actually supports multi-threading at all.
        If ``False``, ``resolve()`` always returns 1 and no flag is emitted.
    """

    name: str = "_thread_spec"
    flag: str | None = "threads"
    max_threads: int | None = None
    fraction: float = 1.0
    multi: bool = True
    dtype: type = str
    default: int = 1
    required = False  # is this required?
    affects_output = False  # do we need a re run?
    description = "_thread_spec for ResourceConfig"  # for documentation
    render = render_value  # "flag", "value", "multi"

    def resolve(self, resources: ResourceConfig) -> int:
        """Compute the concrete thread count for *resources*.

        Parameters
        ----------
        resources : ResourceConfig
            The machine's available resources.

        Returns
        -------
        int
            The thread count to actually use.
        """
        if not self.multi:
            return 1
        n = resources.for_tool(self.fraction)
        if self.max_threads is not None:
            n = min(n, self.max_threads)
        return max(1, n)

    def to_cli(self, resources: ResourceConfig) -> list[str]:
        """Return the CLI fragment for this spec, e.g. ``["-t", "8"]``.

        Returns an empty list when the tool is single-threaded or has no flag.

        Parameters
        ----------
        resources : ResourceConfig
            The machine's available resources.
        """
        if not self.multi or self.flag is None:
            return []
        return [self.flag, str(self.resolve(resources))]


DEFAULT_KEY = "default"


@dataclass(frozen=True)
class ParamSet:
    specs: Mapping[str, ParamSpec]  # key = interner Name, z.B. "threads"
    name: str = DEFAULT_KEY
    subcommand: str | None = None

    def validate(self, params: Params | None) -> None:
        overrides = params.to_dict() if params else {}
        # unknown keys
        unknown = set(overrides) - set(self.specs)
        if unknown:
            raise ValueError(
                f"Unknown params for {self.name} {self.subcommand}: {sorted(unknown)}"
            )

        # required check (default zählt als gesetzt, außer default=None)
        missing: list[str] = []
        for name, spec in self.specs.items():
            value = overrides.get(name, spec.default)
            if spec.required and (value is None or value is False or value == []):
                missing.append(name)
        # if missing:
        #     raise ValueError(
        #         f"Missing required params for {self.name} {self.subcommand}: {missing}"
        #     )  # these are in the arguments

    def merged_values(self, params: Params | None) -> dict[str, Any]:
        overrides = params.to_dict() if params else {}
        merged = {name: spec.default for name, spec in self.specs.items()}
        merged.update(overrides)

        # # normalize
        # for name, spec in self.specs.items():
        #     v = merged.get(name, None)
        #     if v is None:
        #         continue
        #     if spec.normalize:
        #         merged[name] = spec.normalize(v)

        return merged

    def to_cli(self, params: Params | None) -> list[str]:
        if params:
            self.validate(params)
            merged = self.merged_values(params)
            tokens: list[str] = []
            for name, spec in self.specs.items():
                if name in params:
                    v = merged.get(name, None)
                    if v is None or v is False or v == []:
                        continue
                    tokens.extend(spec.render(spec.flag, v))
            return tokens
        return []

    def signature_determinants(self, params: Params | None) -> list[str]:
        """Stable signature tokens for caching/re-run decisions."""
        self.validate(params)
        merged = self.merged_values(params)

        out: list[str] = []
        for name in sorted(self.specs.keys()):
            spec = self.specs[name]
            if not spec.affects_output:
                continue
            v = merged.get(name, None)
            if v is None or v is False or v == []:
                continue
            out.append(f"{name}={v}")
        return out

    def get_spec(self, name: str) -> ParamSpec | None:
        spec = self.specs[name] if name in self.specs else None
        return spec


class ParamRegistry:
    def __init__(
        self, default: ParamSet, by_subcommand: dict[str, ParamSet] | None = None
    ):
        self.default = default
        self.by_subcommand = by_subcommand or {}

    def for_subcommand(self, subcmd: str | None) -> ParamSet:
        if subcmd is None:
            return self.default
        return self.by_subcommand.get(subcmd, self.default)

    def __repr__(self) -> str:
        return (
            f"ParamRegistry(default={self.default}, by_subcommand={self.by_subcommand})"
        )


def initialize_param_registry(
    tool_name: str,
    parameters: Mapping[str, "ParamSet"] | "ParamSet" | str | Path,
) -> ParamRegistry:
    """
    - ParamSet -> default ParamSet
    - Mapping[str, ParamSet] -> uses mapping[DEFAULT_KEY] as default if present
    - str|Path -> JSON file path
    - None -> load resources/params/<tool_name>.json if exists, else empty default
    """
    default = ParamSet({})
    by_subcommand: dict[str, ParamSet] = {}

    if isinstance(parameters, ParamSet):
        return ParamRegistry(default=parameters, by_subcommand={})

    if isinstance(parameters, (str, Path)):
        spec_obj = load_param_json(Path(parameters))
        by_subcommand = _paramsets_from_json(tool_name, spec_obj)
        default = by_subcommand.pop(DEFAULT_KEY, ParamSet({}, tool_name, DEFAULT_KEY))
        return ParamRegistry(default=default, by_subcommand=by_subcommand)

    if isinstance(parameters, Mapping):
        if not all(isinstance(v, ParamSet) for v in parameters.values()):
            raise ValueError("Mapping must be Mapping[str, ParamSet].")
        by_subcommand = dict(parameters)  # shallow copy
        default = by_subcommand.pop(DEFAULT_KEY, ParamSet({}, tool_name, DEFAULT_KEY))
        return ParamRegistry(default=default, by_subcommand=by_subcommand)

    raise ValueError("Invalid parameters format.")


def _thread_spec_from_json(block: dict) -> ToolThreadSpec:
    """Parse a ``_thread_spec`` block from a params JSON file.

    Parameters
    ----------
    block : dict
        The ``_thread_spec`` dict from the JSON, e.g.
        ``{"param_flag": "-t", "fraction": 0.5, "multi": true}``.

    Returns
    -------
    ToolThreadSpec
    """
    return ToolThreadSpec(
        name="_thread_spec",
        flag=block.get("flag", "--threads"),
        max_threads=block.get("max_threads", None),
        fraction=float(block.get("fraction", 1.0)),
        multi=bool(block.get("multi", True)),
    )


def _paramsets_from_json(tool_name: str, obj: dict[str, Any]) -> dict[str, ParamSet]:
    """
    JSON format:
    {
      "mem": {
        "t": {"flag":"-t","dtype":"int","default":1,"required":false,
              "affects_output":false,"render":"value","description":"..."},
        "S": {"flag":"-S","dtype":"bool","default":false,"render":"flag"}
      },
      "index": { ... }
    }
    """
    out: dict[str, ParamSet] = {}

    for subroutine, params_block in obj.items():
        if not isinstance(subroutine, str):
            raise ValueError("Subroutine keys must be strings.")
        if not isinstance(params_block, dict):
            raise ValueError(f"Block for '{subroutine}' must be an object.")

        # Extract _thread_spec before iterating over param specs
        thread_spec: ToolThreadSpec | None = None
        thread_spec_raw = params_block.get("_thread_spec", None)
        if thread_spec_raw is not None:
            thread_spec = _thread_spec_from_json(thread_spec_raw)
            params_block.pop("_thread_spec")
        else:
            thread_spec = ToolThreadSpec(
                name="_thread_spec",
                flag=None,
                max_threads=1,
                fraction=1.0,
                multi=False,
            )

        specs: dict[str, ParamSpec] = {thread_spec.name: thread_spec}

        for key, spec in params_block.items():
            if key == "_thread_spec":  # skip meta key
                continue
            if not isinstance(key, str):
                raise ValueError(f"Param keys in '{subroutine}' must be strings.")
            if not isinstance(spec, dict):
                raise ValueError(f"Spec for '{subroutine}.{key}' must be an object.")

            flag = spec.get("flag", None)
            if flag is not None and not isinstance(flag, str):
                raise ValueError(
                    f"flag must be None or a non-null string value for '{subroutine}.{key}'"
                )

            ps = ParamSpec(
                name=str(spec.get("name", key)),
                flag=flag,
                dtype=_parse_dtype(spec.get("dtype", None)),
                required=bool(spec.get("required", False)),
                default=spec.get("default", None),
                affects_output=bool(spec.get("affects_output", True)),
                description=str(spec.get("description", "")),
                render=_parse_render(spec.get("render", "value")),
                # normalize=_parse_normalize(spec.get("normalize", None)),
            )
            specs[key] = ps

        out[subroutine] = ParamSet(specs, tool_name, subroutine)

    return out


###############################################################################
# Dtypes
###############################################################################

# --- Whitelists for dtypes ---
_DTYPE: dict[str, type] = {
    "str": str,
    "int": int,
    "float": float,
    "bool": bool,
    "any": object,
}


def _parse_dtype(d: Any) -> type | tuple[type, ...]:
    """
    Accepts:
      "int"
      ["int","float"]
    """
    if d is None:
        return object
    if isinstance(d, str):
        if d not in _DTYPE:
            raise ValueError(f"Unknown dtype '{d}'. Allowed: {sorted(_DTYPE)}")
        return _DTYPE[d]
    if isinstance(d, list):
        tys: list[type] = []
        for item in d:
            if not isinstance(item, str) or item not in _DTYPE:
                raise ValueError(
                    f"Unknown dtype item '{item}'. Allowed: {sorted(_DTYPE)}"
                )
            tys.append(_DTYPE[item])
        return tuple(tys)
    raise ValueError(f"Invalid dtype format: {d!r}")


###############################################################################
# Normalize
###############################################################################


# def _parse_normalize(n: Any) -> Callable[[Any], Any] | None:
#     if n is None:
#         return None
#     if isinstance(n, str):
#         if n not in _NORMALIZE:
#             raise ValueError(f"Unknown normalize '{n}'. Allowed: {sorted(_NORMALIZE)}")
#         return _NORMALIZE[n]
#     raise ValueError(f"Invalid normalize format: {n!r}")
