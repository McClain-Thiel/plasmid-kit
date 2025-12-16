from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Dict

import importlib.resources as resources

_CACHE_DIR = Path(os.environ.get("PLASMIDKIT_CACHE", Path.home() / ".cache" / "plasmidkit")).expanduser()
_OFFLINE = bool(int(os.environ.get("PLASMIDKIT_OFFLINE", "0")))


def set_cache_dir(path: str | os.PathLike[str]) -> Path:
    global _CACHE_DIR
    _CACHE_DIR = Path(path).expanduser().absolute()
    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return _CACHE_DIR


def get_cache_dir() -> Path:
    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return _CACHE_DIR


def set_offline(value: bool) -> None:
    global _OFFLINE
    _OFFLINE = bool(value)


def is_offline() -> bool:
    return _OFFLINE


def load_builtin_db(name: str, version: str) -> Dict[str, object]:
    resource_name = "engineered_core_signatures.json"
    if name != "engineered-core" or version != "1.0.0":
        raise FileNotFoundError(f"No built-in database {name}@{version}")
    with resources.files("plasmidkit.data").joinpath(resource_name).open("r", encoding="utf8") as handle:
        return json.load(handle)


def get_artifacts(identifier: str) -> Dict[str, object]:
    if "@" not in identifier:
        raise ValueError("Database identifier must include a version, e.g. name@1.0.0")
    name, version = identifier.split("@", 1)
    return load_builtin_db(name, version)


def ensure_cache_ready() -> None:
    get_cache_dir()
