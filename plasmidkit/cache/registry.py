from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable


@dataclass
class ArtifactSpec:
    url: str
    sha256: str | None = None
    etag: str | None = None


@dataclass
class RegistryEntry:
    name: str
    version: str
    artifacts: Dict[str, ArtifactSpec]


class RegistryManifest:
    def __init__(self, entries: Iterable[RegistryEntry]):
        self._entries = {(entry.name, entry.version): entry for entry in entries}

    def get(self, name: str, version: str) -> RegistryEntry | None:
        return self._entries.get((name, version))

    def list(self) -> Iterable[RegistryEntry]:
        return self._entries.values()
