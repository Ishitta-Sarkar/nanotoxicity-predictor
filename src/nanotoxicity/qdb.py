"""Download and parse the QsarDB zip archives."""

from __future__ import annotations

import hashlib
import urllib.request
import zipfile
from dataclasses import dataclass
from pathlib import Path
from xml.etree import ElementTree

import numpy as np

from nanotoxicity.datasets import FORMULAS, ArchiveSpec

QDB_NS = {"q": "http://www.qsardb.org/QDB"}
USER_AGENT = "nanotoxicity/1.0 (metal-oxide nano-QSAR reanalysis)"


@dataclass(frozen=True)
class OxideTable:
    """Compounds that have both a measured endpoint and every descriptor."""

    key: str
    system: str
    citation: str
    endpoint_label: str
    property_id: str
    descriptor_ids: tuple[str, ...]
    descriptor_names: tuple[str, ...]
    descriptor_units: tuple[str, ...]
    compound_ids: tuple[str, ...]
    names: tuple[str, ...]
    formulas: tuple[str, ...]
    y: np.ndarray
    x: np.ndarray
    excluded_missing_endpoint: tuple[str, ...]


def download_archive(spec: ArchiveSpec, cache_dir: Path) -> Path:
    """Download one archive and reject a file whose SHA-256 does not match."""

    cache_dir.mkdir(parents=True, exist_ok=True)
    destination = cache_dir / spec.filename
    if destination.exists() and _sha256(destination) == spec.sha256:
        return destination
    if destination.exists():
        destination.unlink()

    request = urllib.request.Request(spec.url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=120) as response:
        payload = response.read()
    digest = hashlib.sha256(payload).hexdigest()
    if digest != spec.sha256:
        raise ValueError(
            f"SHA-256 mismatch for {spec.filename}: expected {spec.sha256}, got {digest}"
        )
    destination.write_bytes(payload)
    return destination


def read_archive(path: Path, spec: ArchiveSpec) -> OxideTable:
    """Read descriptors and the endpoint. Compounds without an endpoint are dropped."""

    with zipfile.ZipFile(path) as archive:
        compounds = _compound_names(archive)
        endpoint = _value_table(archive, f"properties/{spec.property_id}/values")
        descriptor_values = []
        descriptor_names = []
        descriptor_units = []
        for descriptor_id in spec.descriptor_ids:
            descriptor_values.append(
                _value_table(archive, f"descriptors/{descriptor_id}/values")
            )
            descriptor_names.append(
                _descriptor_name(archive, f"descriptors/{descriptor_id}")
            )
            descriptor_units.append(
                _unit(archive, f"descriptors/{descriptor_id}/ucum")
            )

    complete_ids = set(compounds)
    for values in descriptor_values:
        complete_ids &= set(values)
    labelled_ids = [cid for cid in endpoint if cid in complete_ids]
    excluded = tuple(sorted(complete_ids - set(endpoint)))

    names = tuple(compounds[cid] for cid in labelled_ids)
    formulas = tuple(_formula(compounds[cid], cid) for cid in labelled_ids)
    y = np.array([endpoint[cid] for cid in labelled_ids], dtype=float)
    columns = [
        np.array([values[cid] for cid in labelled_ids], dtype=float)
        for values in descriptor_values
    ]
    x = np.column_stack(columns) if columns else np.empty((len(labelled_ids), 0))
    return OxideTable(
        key=spec.key,
        system=spec.system,
        citation=spec.citation,
        endpoint_label=spec.endpoint_label,
        property_id=spec.property_id,
        descriptor_ids=spec.descriptor_ids,
        descriptor_names=tuple(descriptor_names),
        descriptor_units=tuple(descriptor_units),
        compound_ids=tuple(labelled_ids),
        names=names,
        formulas=formulas,
        y=y,
        x=x,
        excluded_missing_endpoint=excluded,
    )


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_text(archive: zipfile.ZipFile, name: str) -> str:
    raw = archive.read(name)
    for encoding in ("utf-8-sig", "cp1252"):
        try:
            return raw.decode(encoding)
        except UnicodeDecodeError:
            continue
    return raw.decode("utf-8-sig", errors="replace")


def _compound_names(archive: zipfile.ZipFile) -> dict[str, str]:
    root = ElementTree.fromstring(_read_text(archive, "compounds/compounds.xml"))
    names = {}
    for compound in root.findall("q:Compound", QDB_NS):
        compound_id = compound.findtext("q:Id", default="", namespaces=QDB_NS).strip()
        name = compound.findtext("q:Name", default="", namespaces=QDB_NS).strip()
        names[compound_id] = name
    return names


def _descriptor_name(archive: zipfile.ZipFile, directory: str) -> str:
    root = ElementTree.fromstring(_read_text(archive, "descriptors/descriptors.xml"))
    descriptor_id = directory.rstrip("/").split("/")[-1]
    for descriptor in root.findall("q:Descriptor", QDB_NS):
        if descriptor.findtext("q:Id", default="", namespaces=QDB_NS).strip() == descriptor_id:
            return descriptor.findtext("q:Name", default=descriptor_id, namespaces=QDB_NS).strip()
    return descriptor_id


def _unit(archive: zipfile.ZipFile, name: str) -> str:
    return _read_text(archive, name).strip()


def _value_table(archive: zipfile.ZipFile, name: str) -> dict[str, float]:
    """Parse a QsarDB values file. The first line is a header."""

    values: dict[str, float] = {}
    lines = _read_text(archive, name).splitlines()
    for line in lines[1:]:
        if not line.strip():
            continue
        compound_id, raw_value = line.split("\t", 1)
        raw_value = raw_value.strip()
        if raw_value in {"", "NA", "NaN"}:
            continue
        values[compound_id.strip()] = float(raw_value)
    return values


def _formula(name: str, compound_id: str) -> str:
    if compound_id in FORMULAS:
        return FORMULAS[compound_id]
    return FORMULAS.get(name, compound_id)
