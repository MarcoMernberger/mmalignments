from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterator, Mapping

from mmalignments.services.io import load_param_json

###############################################################################
# Params
###############################################################################


@dataclass(frozen=True)
class Params:
    _override_params: Mapping[str, Any]

    @classmethod
    def of(cls, **kwargs: Any) -> "Params":
        return cls(kwargs)

    def __init__(self, **kwargs: Any) -> None:
        object.__setattr__(self, "_override_params", dict(kwargs))

    # ---- Mapping interface ----
    def __getitem__(self, key: str) -> Any:
        return self._override_paramss[key]

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

    def add_new(self, **kwargs: Any) -> "Params":
        return Params(**self.to_dict(), **kwargs)

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
    required: True | False = False  # is this required?
    default: Any = None  # default value
    affects_output: bool = True  # do we need a re run?
    description: str = ""  # for documentation
    render: Callable = render_value  # "flag", "value", "multi"


#    normalize: Callable[[Any], Any] | None = None  # optional normalization function

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
    parameters: Mapping[str, "ParamSet"] | "ParamSet" | str | Path | None,
) -> ParamRegistry:
    """
    - ParamSet -> default ParamSet
    - Mapping[str, ParamSet] -> uses mapping[DEFAULT_KEY] as default if present
    - str|Path -> JSON file path
    - None -> load resources/params/<tool_name>.json if exists, else empty default
    """
    default = ParamSet({})
    by_subcommand: dict[str, ParamSet] = {}

    if parameters is None:
        # try implicit default json
        default_path = (
            Path(__file__).resolve().parent
            / "resources"
            / "params"
            / f"{tool_name}.json"
        )
        if default_path.exists():
            spec_obj = load_param_json(default_path)
            by_subcommand = _paramsets_from_json(tool_name, spec_obj)
            default = by_subcommand.pop(
                DEFAULT_KEY, ParamSet({}, tool_name, DEFAULT_KEY)
            )
        return ParamRegistry(default=default, by_subcommand=by_subcommand)

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

        specs: dict[str, ParamSpec] = {}

        for key, spec in params_block.items():
            if not isinstance(key, str):
                raise ValueError(f"Param keys in '{subroutine}' must be strings.")
            if not isinstance(spec, dict):
                raise ValueError(f"Spec for '{subroutine}.{key}' must be an object.")

            flag = spec.get("flag", None)
            if flag is not None and not isinstance(flag, str):
                raise ValueError(
                    f"flag must be string or null for '{subroutine}.{key}'"
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
