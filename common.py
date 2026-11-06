"""Shared helpers: download + cache GenBank records from NCBI, and a common plot style."""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # write PNG files; no window needed
import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO

BASE = Path(__file__).resolve().parent
DATA = BASE / "data"
RESULTS = BASE / "results"

INK = "#2B2D31"
ACCENT = "#B5562F"
SOFT = "#8A8F98"


def fetch_genbank(accession: str):
    """Return a SeqRecord for *accession*; the file is cached in ./data after the first download."""
    DATA.mkdir(exist_ok=True)
    path = DATA / f"{accession}.gb"
    if not path.exists():
        email = os.environ.get("ENTREZ_EMAIL", "").strip()
        if not email:
            raise SystemExit(
                "NCBI asks for a contact e-mail. Set it first, e.g. on Windows:\n"
                "    set ENTREZ_EMAIL=you@example.com\n"
                "then run the script again."
            )
        Entrez.email = email
        print(f"Downloading {accession} from NCBI ...")
        with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as handle:
            path.write_text(handle.read(), encoding="utf-8", newline="\n")
    return SeqIO.read(path, "genbank")


def cds_feature(record, index: int = 0):
    """Return the *index*-th CDS feature of a record."""
    cds = [f for f in record.features if f.type == "CDS"]
    if not cds:
        raise SystemExit(f"{record.id} has no CDS feature.")
    return cds[index]


def aa_percent(analysis) -> dict:
    """Amino-acid percentages (0-100) from a ProteinAnalysis, for old and new Biopython."""
    value = getattr(analysis, "amino_acids_percent", None)
    if value is None:                       # Biopython < 1.80
        value = analysis.get_amino_acids_percent()
    elif callable(value):
        value = value()
    return {aa: v * 100 for aa, v in value.items()}


def set_style() -> None:
    plt.rcParams.update({
        "figure.dpi": 120,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": "#E3E3E3",
        "axes.axisbelow": True,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "text.color": INK,
        "font.size": 10,
        "axes.titleweight": "bold",
    })


def save(fig, name: str) -> None:
    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / name
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"saved {out.relative_to(BASE)}")
