from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict

import importlib.resources as resources

_CACHE_DIR = Path(os.environ.get("PLASMIDKIT_CACHE", Path.home() / ".cache" / "plasmidkit")).expanduser()
_OFFLINE = bool(int(os.environ.get("PLASMIDKIT_OFFLINE", "0")))

_REGISTRIES: Dict[str, "Registry"] = {}


@dataclass
class Artifact:
    name: str
    path: Path
    sha256: str | None = None


class Registry:
    def __init__(self, name: str, manifest_path: Path):
        self.name = name
        self.manifest_path = manifest_path
        self._manifest_data: Dict[str, object] | None = None

    def load(self) -> Dict[str, object]:
        if self._manifest_data is None:
            with open(self.manifest_path, "r", encoding="utf8") as handle:
                import yaml  # local import to avoid dependency if unused

                self._manifest_data = yaml.safe_load(handle)
        return self._manifest_data


def add_registry(name: str, manifest: str | os.PathLike[str]) -> Registry:
    path = Path(manifest).expanduser()
    if not path.exists():
        raise FileNotFoundError(path)
    registry = Registry(name, path)
    _REGISTRIES[name] = registry
    return registry


def list_registries() -> Dict[str, Registry]:
    return dict(_REGISTRIES)


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
    try:
        data = load_builtin_db(name, version)
    except FileNotFoundError as exc:  # pragma: no cover - placeholder for future registry support
        raise RuntimeError(str(exc)) from exc
    return data


def ensure_cache_ready() -> None:
    get_cache_dir()
