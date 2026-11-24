import hashlib
import io
import zipfile
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from nanotoxicity.datasets import GAJEWICZ
from nanotoxicity.qdb import download_archive, read_archive


def _archive_bytes() -> bytes:
    compounds = """<?xml version="1.0" encoding="UTF-8"?>
<CompoundRegistry xmlns="http://www.qsardb.org/QDB">
  <Compound><Id>ZnO</Id><Name>Zinc II oxide</Name></Compound>
  <Compound><Id>TiO2</Id><Name>Titanium IV oxide</Name></Compound>
  <Compound><Id>P_1</Id><Name>Copper II oxide</Name></Compound>
</CompoundRegistry>
"""
    descriptors = """<?xml version="1.0" encoding="UTF-8"?>
<DescriptorRegistry xmlns="http://www.qsardb.org/QDB">
  <Descriptor><Id>dH-fc</Id><Name>Enthalpy</Name></Descriptor>
  <Descriptor><Id>Xc</Id><Name>Electronegativity</Name></Descriptor>
</DescriptorRegistry>
"""
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("compounds/compounds.xml", compounds)
        archive.writestr("descriptors/descriptors.xml", descriptors)
        archive.writestr(
            "descriptors/dH-fc/values",
            "Compound Id\tEnthalpy\nZnO\t-10\nTiO2\t-20\nP_1\t-30\n",
        )
        archive.writestr("descriptors/dH-fc/ucum", "kcal/mol\n")
        archive.writestr(
            "descriptors/Xc/values",
            "Compound Id\tElectronegativity\nZnO\t1.5\nTiO2\t2.5\nP_1\t3.5\n",
        )
        archive.writestr("descriptors/Xc/ucum", "eV\n")
        archive.writestr(
            "properties/pLC50/values",
            "Compound Id\tpLC50\nZnO\t3.0\nTiO2\t2.0\n",
        )
        archive.writestr(
            "predictions/pLC50_train/values",
            "Compound Id\tprediction\nZnO\t99\nTiO2\t99\n",
        )
    return buffer.getvalue()


def test_cp1252_descriptor_header_is_readable(tmp_path: Path):
    # The Gajewicz electronegativity header uses a Windows-1252 apostrophe (byte 0x92).
    blob = _archive_bytes()
    with zipfile.ZipFile(io.BytesIO(blob)) as original, zipfile.ZipFile(
        tmp_path / "cp1252.zip", "w"
    ) as rewritten:
        for info in original.infolist():
            data = original.read(info.filename)
            if info.filename == "descriptors/Xc/values":
                data = data.replace(b"Electronegativity", b"Mulliken\x92s electronegativity")
            rewritten.writestr(info.filename, data)
    table = read_archive(tmp_path / "cp1252.zip", GAJEWICZ)
    assert table.descriptor_names[1] == "Electronegativity"
    np.testing.assert_allclose(table.x[:, 1], [1.5, 2.5])
    path = tmp_path / "toy.zip"
    path.write_bytes(_archive_bytes())
    table = read_archive(path, GAJEWICZ)
    assert table.compound_ids == ("ZnO", "TiO2")
    assert table.formulas == ("ZnO", "TiO2")
    assert table.excluded_missing_endpoint == ("P_1",)
    assert table.descriptor_names == ("Enthalpy", "Electronegativity")
    assert table.descriptor_units == ("kcal/mol", "eV")
    np.testing.assert_allclose(table.y, [3.0, 2.0])
    np.testing.assert_allclose(table.x, [[-10.0, 1.5], [-20.0, 2.5]])
    assert not np.allclose(table.x[:, 0], table.y)


def test_checksum_mismatch_is_rejected(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    payload = b"not-the-archive"
    monkeypatch.setattr(
        "nanotoxicity.qdb.urllib.request.urlopen",
        lambda request, timeout=120: _Response(payload),
    )
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        download_archive(GAJEWICZ, tmp_path)


def test_matching_checksum_is_cached(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    payload = b"cached-bytes"
    spec = replace(GAJEWICZ, sha256=hashlib.sha256(payload).hexdigest())
    calls = {"count": 0}

    def fake_open(request, timeout=120):
        calls["count"] += 1
        return _Response(payload)

    monkeypatch.setattr("nanotoxicity.qdb.urllib.request.urlopen", fake_open)
    first = download_archive(spec, tmp_path)
    second = download_archive(spec, tmp_path)
    assert first == second
    assert calls["count"] == 1
    assert first.read_bytes() == payload


class _Response:
    def __init__(self, payload: bytes):
        self._payload = payload

    def read(self) -> bytes:
        return self._payload

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        return False
