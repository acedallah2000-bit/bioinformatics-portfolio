"""Project 3 - Comparing the genes and proteins of SARS-CoV-2 (reference genome NC_045512.2).

What it does
  * extracts every CDS from the GenBank record and translates it
  * bar charts: protein length and GC content per gene
  * heat map: amino-acid composition of each protein
  * Kyte-Doolittle hydropathy of the Spike protein (signal peptide / transmembrane region)
  * table with length, GC%, molecular weight, pI and GRAVY per protein
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.SeqUtils import ProtParamData, gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from common import ACCENT, INK, RESULTS, SOFT, aa_percent, fetch_genbank, save, set_style

ACCESSION = "NC_045512.2"
AA_ORDER = "AVLIMFWCGPYSTNQHKRDE"   # hydrophobic -> special -> polar -> charged


def main() -> None:
    set_style()
    record = fetch_genbank(ACCESSION)
    genome_gc = gc_fraction(str(record.seq)) * 100
    print(f"{record.id}: {len(record.seq):,} nt, GC {genome_gc:.1f}%")

    rows, composition, proteins = [], {}, {}
    for feat in record.features:
        if feat.type != "CDS":
            continue
        gene = feat.qualifiers.get("gene", ["?"])[0]
        product = feat.qualifiers.get("product", [""])[0]
        if gene == "ORF1ab" and "ORF1ab" not in product:
            continue  # ORF1a is nested inside ORF1ab - keep only the longer one
        nt = feat.extract(record.seq)
        prot = str(nt.translate(to_stop=True))
        pa = ProteinAnalysis(prot)
        proteins[gene] = prot
        composition[gene] = aa_percent(pa)
        rows.append({
            "gene": gene, "product": product, "length_nt": len(nt), "length_aa": len(prot),
            "gc_percent": round(gc_fraction(str(nt)) * 100, 2),
            "mw_kda": round(pa.molecular_weight() / 1000, 2),
            "isoelectric_point": round(pa.isoelectric_point(), 2),
            "gravy": round(pa.gravy(), 3),
        })
    table = pd.DataFrame(rows)
    RESULTS.mkdir(exist_ok=True)
    table.to_csv(RESULTS / "03_sars_cov2_proteins.csv", index=False)
    print(table[["gene", "length_aa", "gc_percent", "mw_kda"]].to_string(index=False))

    # 1) length + GC per gene
    by_len = table.sort_values("length_aa")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4))
    a1.barh(by_len["gene"], by_len["length_aa"], color=INK)
    a1.set_xscale("log")
    a1.set_xlabel("Protein length (aa, log scale)")
    a1.set_title("Protein length")
    a1.grid(axis="y", visible=False)
    a2.barh(table["gene"][::-1], table["gc_percent"][::-1], color=ACCENT)
    a2.axvline(genome_gc, color=SOFT, linestyle=":", label=f"genome mean {genome_gc:.1f}%")
    a2.set_xlabel("GC %")
    a2.set_title("GC content per gene")
    a2.legend(frameon=False, loc="lower right", fontsize=8)
    a2.grid(axis="y", visible=False)
    fig.tight_layout()
    save(fig, "03_length_and_gc.png")

    # 2) amino-acid heat map
    genes = list(table["gene"])
    matrix = np.array([[composition[g].get(aa, 0.0) for aa in AA_ORDER] for g in genes])
    fig, ax = plt.subplots(figsize=(9, 4.6))
    image = ax.imshow(matrix, aspect="auto", cmap="Greys")
    ax.set_xticks(range(len(AA_ORDER)), labels=list(AA_ORDER))
    ax.set_yticks(range(len(genes)), labels=genes)
    ax.grid(False)
    ax.set_title("Amino-acid composition of SARS-CoV-2 proteins")
    fig.colorbar(image, ax=ax, label="% of residues")
    save(fig, "03_amino_acid_heatmap.png")

    # 3) Spike hydropathy
    spike = proteins["S"]
    window = 19
    scale = ProteinAnalysis(spike).protein_scale(ProtParamData.kd, window=window)
    pos = [i + (window + 1) / 2 for i in range(len(scale))]
    fig, ax = plt.subplots(figsize=(9, 3.2))
    ax.plot(pos, scale, color=INK, linewidth=1.1)
    ax.fill_between(pos, scale, 1.6, where=[s > 1.6 for s in scale], color=ACCENT, alpha=0.5,
                    label="> 1.6 (possible membrane-spanning segment)")
    ax.axhline(0, color=INK, linewidth=0.8)
    ax.set_xlabel("Residue position")
    ax.set_ylabel(f"Hydropathy (window {window})")
    ax.set_title(f"Spike protein hydropathy ({len(spike)} aa)")
    ax.legend(frameon=False, loc="lower left", fontsize=8)
    save(fig, "03_spike_hydropathy.png")


if __name__ == "__main__":
    main()
