"""Project 1 - Anatomy of a human gene: HBB (beta-globin) mRNA, NCBI NM_000518.

What it does
  * downloads the mRNA record and extracts the coding sequence (CDS)
  * GC content along the whole mRNA (sliding window), with the CDS highlighted
  * translates the CDS and plots the amino-acid composition
  * Kyte-Doolittle hydropathy profile of the protein
  * writes a one-row summary table (length, GC%, molecular weight, pI, GRAVY)
"""
import pandas as pd
import matplotlib.pyplot as plt
from Bio.SeqUtils import ProtParamData, gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from common import (ACCENT, INK, RESULTS, SOFT, aa_percent, cds_feature, fetch_genbank, save,
                    set_style)

ACCESSION = "NM_000518"


def sliding_gc(seq: str, window: int = 50, step: int = 5):
    """GC % of overlapping windows; x = window centre (nt)."""
    xs, ys = [], []
    for i in range(0, len(seq) - window + 1, step):
        xs.append(i + window / 2)
        ys.append(gc_fraction(seq[i:i + window]) * 100)
    return xs, ys


def main() -> None:
    set_style()
    record = fetch_genbank(ACCESSION)
    mrna = str(record.seq).upper()
    feat = cds_feature(record)
    cds = feat.extract(record.seq)
    protein = str(cds.translate(to_stop=True))
    analysis = ProteinAnalysis(protein)
    print(f"{record.id}: {record.description}")
    print(f"mRNA {len(mrna)} nt | CDS {len(cds)} nt | protein {len(protein)} aa")

    # 1) GC content along the mRNA
    xs, ys = sliding_gc(mrna)
    mean_gc = gc_fraction(mrna) * 100
    fig, ax = plt.subplots(figsize=(8, 3.2))
    ax.axvspan(int(feat.location.start), int(feat.location.end), color=ACCENT, alpha=0.12,
               label="coding sequence (CDS)")
    ax.plot(xs, ys, color=INK, linewidth=1.3)
    ax.axhline(mean_gc, color=SOFT, linestyle=":", label=f"mean {mean_gc:.1f}%")
    ax.set_xlabel("Position in mRNA (nt)")
    ax.set_ylabel("GC % (50-nt window)")
    ax.set_title(f"GC content of {record.id}")
    ax.legend(frameon=False, loc="lower right")
    save(fig, "01_hbb_gc_content.png")

    # 2) amino-acid composition
    pct = aa_percent(analysis)
    ordered = sorted(pct.items(), key=lambda kv: kv[1], reverse=True)
    fig, ax = plt.subplots(figsize=(8, 3.2))
    ax.bar([a for a, _ in ordered], [v for _, v in ordered], color=INK, width=0.7)
    ax.set_ylabel("Frequency (%)")
    ax.set_title(f"Amino-acid composition of beta-globin ({len(protein)} aa)")
    save(fig, "01_hbb_amino_acids.png")

    # 3) Kyte-Doolittle hydropathy
    window = 9
    scale = analysis.protein_scale(ProtParamData.kd, window=window)
    pos = [i + (window + 1) / 2 for i in range(len(scale))]
    fig, ax = plt.subplots(figsize=(8, 3.2))
    ax.plot(pos, scale, color=INK, linewidth=1.3)
    ax.fill_between(pos, scale, 0, where=[s > 0 for s in scale], color=ACCENT, alpha=0.3,
                    label="hydrophobic")
    ax.fill_between(pos, scale, 0, where=[s <= 0 for s in scale], color=SOFT, alpha=0.3,
                    label="hydrophilic")
    ax.axhline(0, color=INK, linewidth=0.8)
    ax.set_xlabel("Residue position")
    ax.set_ylabel(f"Hydropathy (window {window})")
    ax.set_title("Kyte-Doolittle hydropathy profile")
    ax.legend(frameon=False)
    save(fig, "01_hbb_hydropathy.png")

    # 4) summary table
    summary = pd.DataFrame([{
        "accession": record.id,
        "mrna_length_nt": len(mrna),
        "cds_length_nt": len(cds),
        "protein_length_aa": len(protein),
        "mrna_gc_percent": round(mean_gc, 2),
        "cds_gc_percent": round(gc_fraction(str(cds)) * 100, 2),
        "protein_mw_da": round(analysis.molecular_weight(), 1),
        "isoelectric_point": round(analysis.isoelectric_point(), 2),
        "gravy": round(analysis.gravy(), 3),
    }])
    summary.to_csv(RESULTS / "01_hbb_summary.csv", index=False)
    print(summary.T.to_string(header=False))


if __name__ == "__main__":
    main()
