import os
from pathlib import Path

def get_variable(var_name: str) -> str:
    return os.environ[var_name]

def hostname() -> str:
    return os.environ["NICE_HOSTNAME"]

def prebuilt_path() -> Path:
    return Path(os.environ["MBF_EXTERNAL_PREBUILD_PATH"])