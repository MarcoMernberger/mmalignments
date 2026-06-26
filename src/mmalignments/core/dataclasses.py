from typing import Any


class PrefixDict(dict):

    def __init__(self, data=None, separator="_", **kwargs):
        if data is None:
            data = {}
        self.separator = separator
        super().__init__({**data, **kwargs})
        self._validate()

    def _validate(self):
        keys = set(self.keys())

        for k in keys:
            # case 1: k is namespace root
            for other in keys:
                if other.startswith(k + self.separator) and other != k:
                    raise ValueError(
                        f"Ambiguous namespace: '{k}' is prefix of other keys: {other}"
                    )

                if other.endswith(self.separator + k) and other != k:
                    raise ValueError(
                        f"Ambiguous namespace: '{k}' is suffix of other keys: {other}"
                    )

    def group(self, prefix: str) -> dict[str, Any]:
        """Explicitly group keys by prefixes."""
        return self.prefix(prefix)

    def prefix(self, prefix: str) -> dict[str, Any]:
        """Get a sub-dictionary of all keys that start with the given prefix."""
        prefix_key = f"{prefix}{self.separator}"

        return {k: v for k, v in self.items() if k.startswith(prefix_key)}

    def __getattr__(self, key: str) -> Any:
        """Access keys as attributes, with support for prefix namespaces."""
        # 1. exact key wins
        if key in self:
            return self[key]

        # 2. fallback: prefix namespace
        sub = self.prefix(key)
        if sub:
            return sub

        raise AttributeError(f"No key or prefix namespace '{key}'")

    def update(self, *args, **kwargs):
        super().update(*args, **kwargs)
        self._validate()
