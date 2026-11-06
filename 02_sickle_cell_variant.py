"""Project 2 - From one DNA letter to one amino acid: the sickle-cell variant of HBB.

The reference beta-globin CDS (NM_000518) is edited in silico at c.20A>T, the well-known
HbS variant (rs334). The script then
  * translates both versions and reports the protein change
  * compares hydropathy profiles (wild type vs variant)
  * compares isoelectric point, GRAVY and molecular weight
"""
import pandas as pd
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.SeqUtils import ProtParamData, seq3
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from common import ACCENT, INK, RESULTS, SOFT, cds_feature, fetch_genbank, save, set_style

ACCESSION = "NM_000518"
POSITION = 19            # 0-based index of c.20 inside the CDS
EXPECTED_CODON = "GAG"   # Glu


def main() -> None:
    set_style()
    record = fetch_genbank(ACCESSION)
    cds = str(cds_feature(record).extract(record.seq)).upper()

    # Safety check: make sure we are editing the right base.
    if cds[POSITION - 1:POSITION + 2] != EXPECTED_CODON:
        raise SystemExit(f"Unexpected codon {cds[POSITION - 1:POSITION + 2]} at c.20 - check the record.")
    variant = cds[:POSITION] + "T" + cds[POSITION + 1:]

    wt_prot = str(Seq(cds).translate(to_stop=True))
    mut_prot = str(Seq(variant).translate(to_stop=True))
    changes = [(i + 1, a, b) for i, (a, b) in enumerate(zip(wt_prot, mut_prot)) if a != b]
    print(f"DNA change : c.20A>T  ({cds[POSITION - 1:POSITION + 2]} -> {variant[POSITION - 1:POSITION + 2]})")
    for pos, a, b in changes:
        print(f"Protein    : p.{seq3(a)}{pos}{seq3(b)}  (legacy numbering without Met: {seq3(a)}{pos - 1}{seq3(b)})")
    print(f"Length     : {len(wt_prot)} aa (wild type) vs {len(mut_prot)} aa (variant)")

    wt, mut = ProteinAnalysis(wt_prot), ProteinAnalysis(mut_prot)

    # 1) hydropathy overlay
    window = 7
    s_wt = wt.protein_scale(ProtParamData.kd, window=window)
    s_mut = mut.protein_scale(ProtParamData.kd, window=window)
    pos = [i + (window + 1) / 2 for i in range(len(s_wt))]
    fig, ax = plt.subplots(figsize=(8, 3.4))
    ax.plot(pos, s_wt, color=INK, linewidth=1.4, label="wild type (HbA)")
    ax.plot(pos, s_mut, color=ACCENT, linewidth=1.4, linestyle="--", label="variant (HbS)")
    for p, _, _ in changes:
        ax.axvline(p, color=SOFT, linestyle=":", linewidth=1)
        ax.annotate(f"residue {p}", (p, max(s_wt + s_mut)), xytext=(4, -2),
                    textcoords="offset points", fontsize=8, color=SOFT)
    ax.axhline(0, color=INK, linewidth=0.8)
    ax.set_xlabel("Residue position")
    ax.set_ylabel(f"Hydropathy (window {window})")
    ax.set_title("Hydropathy: wild type vs sickle-cell variant")
    ax.legend(frameon=False, loc="lower right")
    save(fig, "02_hydropathy_wt_vs_sickle.png")

    # 2) property comparison
    props = {
        "Isoelectric point (pI)": (wt.isoelectric_point(), mut.isoelectric_point()),
        "GRAVY (hydropathy)": (wt.gravy(), mut.gravy()),
        "Molecular weight (kDa)": (wt.molecular_weight() / 1000, mut.molecular_weight() / 1000),
    }
    fig, axes = plt.subplots(1, 3, figsize=(9, 3))
    for ax, (name, (a, b)) in zip(axes, props.items()):
        ax.bar(["WT", "HbS"], [a, b], color=[INK, ACCENT], width=0.55)
        ax.set_title(name, fontsize=9)
        ax.grid(axis="x", visible=False)
        for x, v in enumerate((a, b)):
            ax.text(x, v, f"{v:.3f}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    save(fig, "02_property_changes.png")

    # 3) variant table
    rows = [{
        "dna_change": "c.20A>T",
        "codon_change": f"{EXPECTED_CODON}>{variant[POSITION - 1:POSITION + 2]}",
        "protein_change": f"p.{seq3(a)}{p}{seq3(b)}",
        "wt_pI": round(wt.isoelectric_point(), 3), "variant_pI": round(mut.isoelectric_point(), 3),
        "wt_gravy": round(wt.gravy(), 3), "variant_gravy": round(mut.gravy(), 3),
        "wt_mw_da": round(wt.molecular_weight(), 1), "variant_mw_da": round(mut.molecular_weight(), 1),
    } for p, a, b in changes]
    pd.DataFrame(rows).to_csv(RESULTS / "02_variant_summary.csv", index=False)


if __name__ == "__main__":
    main()
