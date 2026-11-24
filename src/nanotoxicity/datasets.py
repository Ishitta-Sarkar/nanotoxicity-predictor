"""QsarDB archives read by this reanalysis.

The checksums are SHA-256 of the zip bytes. They pin the files that were
downloaded while preparing the reported results.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ArchiveSpec:
    """One published QsarDB zip and the columns used here."""

    key: str
    filename: str
    url: str
    sha256: str
    property_id: str
    descriptor_ids: tuple[str, ...]
    endpoint_label: str
    system: str
    citation: str
    doi: str
    data_doi: str


GAJEWICZ = ArchiveSpec(
    key="gajewicz_hacat",
    filename="2015N313.qdb.zip",
    url="https://qsardb.org/repository/bitstream/handle/10967/214/2015N313.qdb.zip",
    sha256="8888e614436e8a77faecb2d5d72624e8dc00f5a6023c2e699537679ee36fec77",
    property_id="pLC50",
    descriptor_ids=("dH-fc", "Xc"),
    endpoint_label="log(1/LC50), with LC50 in mol/L",
    system="HaCaT human keratinocytes",
    citation="Gajewicz et al., 2014",
    doi="10.3109/17435390.2014.930195",
    data_doi="10.15152/QDB.214",
)

PUZYN = ArchiveSpec(
    key="puzyn_ecoli",
    filename="2011NN175.zip",
    url="https://qsardb.org/repository/bitstream/handle/10967/119/2011NN175.zip",
    sha256="7e16d28c40849b6937aab3221b3a726891a39aefb9b33bb9ef60291a189ccfa2",
    property_id="pEC50",
    descriptor_ids=("dH_Me",),
    endpoint_label="log(1/EC50), with EC50 in mol/L",
    system="Escherichia coli",
    citation="Puzyn et al., 2011",
    doi="10.1038/nnano.2011.10",
    data_doi="10.15152/QDB.119",
)

ARCHIVES: tuple[ArchiveSpec, ...] = (GAJEWICZ, PUZYN)

# Display formulas for the compound names stored in the Gajewicz archive.
# Puzyn compound ids are already formulas.
FORMULAS = {
    "Titanium IV oxide": "TiO2",
    "Aluminum III oxide": "Al2O3",
    "Zirconium IV oxide": "ZrO2",
    "Iron III oxide": "Fe2O3",
    "Silicon IV dioxide": "SiO2",
    "Yttrium III oxide": "Y2O3",
    "Vanadium III oxide": "V2O3",
    "Chromium III oxide": "Cr2O3",
    "Antimony III oxide": "Sb2O3",
    "Nickel II oxide": "NiO",
    "Bismuth III oxide": "Bi2O3",
    "Tungsten VI oxide": "WO3",
    "Manganese III oxide": "Mn2O3",
    "Tin IV oxide": "SnO2",
    "Cobalt II oxide": "CoO",
    "Lanthanum III oxide": "La2O3",
    "Indium III oxide": "In2O3",
    "Zinc II oxide": "ZnO",
    "Copper II oxide": "CuO",
    "Iron II oxide": "FeO",
    "Gadolinium III oxide": "Gd2O3",
    "Lead II oxide": "PbO",
    "Lead IV oxide": "PbO2",
}
