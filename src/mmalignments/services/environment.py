import os
from pathlib import Path


def get_variable(var_name: str) -> str:
    return os.environ[var_name]


def hostname() -> str:
    if "NICE_HOSTNAME" not in os.environ:
        return "no nice hostbname set"
    return os.environ["NICE_HOSTNAME"]


def prebuilt_path() -> Path:
    if "MBF_EXTERNAL_PREBUILD_PATH" not in os.environ:
        return ""
    return Path(os.environ["MBF_EXTERNAL_PREBUILD_PATH"])
