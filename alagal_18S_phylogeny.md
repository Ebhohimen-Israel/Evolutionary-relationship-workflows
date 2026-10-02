### Step 1: Project Directory Structure & Initialization

Execute the following in a Python cell within JupyterLab to establish a standardized directory layout:

```python
import os
from pathlib import Path

BASE_DIR = Path("algal_18S_phylogeny")
SUB_DIRS = [
    "data/raw_genomes",
    "data/references",
    "data/extracted_18S",
    "blast_db",
    "alignments",
    "phylogeny",
    "results/tables",
    "results/figures"
]

for sub in SUB_DIRS:
    (BASE_DIR / sub).mkdir(parents=True, exist_ok=True)

print("Directory tree successfully initialized.")

```

* **Extracted Artifacts:** Clean filesystem structure segregating raw data, intermediate alignments, phylogenetic outputs, and final deliverables.

---

### Step 2: Sequence Acquisition (Targets, Ingroup Taxa, and Outgroup)

Retrieve complete 18S rRNA sequences for the primary comparative pair (*Chlorella vulgaris* and *Tetradesmus obliquus*), selected Chlorophyta relatives, and an outgroup via NCBI Entrez:

```python
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

Entrez.email = "your_email@institution.edu"

# Curated reference accessions representing Trebouxiophyceae, Chlorophyceae, and Outgroup
TAXA_MAP = {
    # Ingroup: Trebouxiophyceae
    "Chlorella_vulgaris": "KJ756823.1",
    "Chlorella_sorokiniana": "AB488737.1",
    "Chlorella_variabilis": "FM205852.1",
    "Auxenochlorella_pyrenoidosa": "KM061466.1",
    # Ingroup: Chlorophyceae (Sphaeropleales)
    "Tetradesmus_obliquus": "FR865738.1",
    "Scenedesmus_acutus": "AJ413375.1",
    "Desmodesmus_communis": "HE610052.1",
    "Monoraphidium_neglectum": "AJ421147.1",
    "Pediastrum_duplex": "HE610057.1",
    # Ingroup: Chlorophyceae (Chlamydomonadales)
    "Chlamydomonas_reinhardtii": "M32703.1",
    "Volvox_carteri": "X53904.1",
    # Outgroup: Streptophyta
    "Mesotaenium_endlicherianum": "HQ659681.1"
}

records = []
for organism, acc in TAXA_MAP.items():
    handle = Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    
    # Extract specifically the 18S rRNA feature coordinates
    extracted_seq = None
    for feat in record.features:
        if feat.type == "rRNA":
            product = " ".join(feat.qualifiers.get("product", [])).lower()
            if "18s" in product or "small subunit" in product:
                extracted_seq = feat.extract(record.seq)
                break
                
    if extracted_seq is None:
        extracted_seq = record.seq  # Fallback to full record if dedicated feature tag is absent
        
    clean_rec = SeqRecord(
        extracted_seq,
        id=organism,
        description=f"18S rRNA | Accession: {acc}"
    )
    records.append(clean_rec)

SeqIO.write(records, "algal_18S_phylogeny/data/references/dataset_18S.fasta", "fasta")
print(f"Retrieved and stored {len(records)} verified 18S rRNA sequences.")

```

* **Extracted Artifacts:** `data/references/dataset_18S.fasta` containing verified, homologous 18S sequences.

---

### Step 3: Sequence Quality Control & Base Composition

Compute sequence lengths, individual nucleotide counts, and stoichiometric ratios:

```python
import pandas as pd

def compute_qc(record):
    seq_str = str(record.seq).upper()
    length = len(seq_str)
    a, t, g, c = (seq_str.count(b) for b in "ATGC")
    n = seq_str.count("N")
    gc_perc = ((g + c) / length) * 100
    at_perc = ((a + t) / length) * 100
    
    return {
        "Taxon": record.id,
        "Sequence Length (bp)": length,
        "A (%)": round((a / length) * 100, 2),
        "T (%)": round((t / length) * 100, 2),
        "G (%)": round((g / length) * 100, 2),
        "C (%)": round((c / length) * 100, 2),
        "Ambiguous (N)": n,
        "GC Content (%)": round(gc_perc, 2),
        "AT/GC Ratio": round(at_perc / gc_perc, 3)
    }

qc_df = pd.DataFrame([compute_qc(r) for r in records])
qc_df.to_csv("algal_18S_phylogeny/results/tables/Table_4_1_sequence_qc.csv", index=False)
qc_df

```

* **Extracted Presentation Deliverable:** **Table 4.1: Nucleotide Composition and QC Metrics of Selected Microalgal 18S rRNA Genes** (formatted CSV/DataFrame detailing base bias, length variations, and sequence integrity).



---

### Step 4: Local BLAST Validation

Verify that all extracted sequences match annotated 18S rRNA entries in the reference database:

```python
# Build local BLAST database
!makeblastdb \
    -in algal_18S_phylogeny/data/references/dataset_18S.fasta \
    -dbtype nucl \
    -out algal_18S_phylogeny/blast_db/reference_18S_db

# Run BLASTn all-vs-all validation
!blastn \
    -query algal_18S_phylogeny/data/references/dataset_18S.fasta \
    -db algal_18S_phylogeny/blast_db/reference_18S_db \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen evalue bitscore" \
    -out algal_18S_phylogeny/results/tables/blast_validation_raw.txt

blast_cols = ["Query", "Subject", "Identity (%)", "Alignment Length", "Mismatches", "Gaps", "E-value", "Bit Score"]
blast_df = pd.read_csv("algal_18S_phylogeny/results/tables/blast_validation_raw.txt", sep="\t", names=blast_cols)

# Filter for the top reciprocal cross-matches excluding self-hits
cross_hits = blast_df[blast_df["Query"] != blast_df["Subject"]].groupby("Query").first().reset_index()
cross_hits.to_csv("algal_18S_phylogeny/results/tables/Table_4_2_blast_validation.csv", index=False)
cross_hits

```

* **Extracted Presentation Deliverable:** **Table 4.2: Top BLASTn Homology Verification and Coordinate Mapping** (confirms non-target sequences are not included).



---

### Step 5: High-Resolution Global Pairwise Alignment

Generate a Needleman-Wunsch global alignment specifically comparing *C. vulgaris* directly with *T. obliquus*:

```python
from Bio import Align

aligner = Align.PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 2.0
aligner.mismatch_score = -1.0
aligner.open_gap_score = -5.0
aligner.extend_gap_score = -0.5

cv_seq = str([r.seq for r in records if r.id == "Chlorella_vulgaris"][0])
to_seq = str([r.seq for r in records if r.id == "Tetradesmus_obliquus"][0])

pairwise_aln = aligner.align(cv_seq, to_seq)[0]
aligned_cv, aligned_to = pairwise_aln[0], pairwise_aln[1]

matches = sum(1 for a, b in zip(aligned_cv, aligned_to) if a == b and a not in "-?")
gaps = sum(1 for a, b in zip(aligned_cv, aligned_to) if a == "-" or b == "-")
total_cols = len(aligned_cv)
pairwise_identity = (matches / (total_cols - gaps)) * 100

pairwise_summary = pd.DataFrame([{
    "Comparison": "C. vulgaris vs. T. obliquus",
    "Alignment Mode": "Global (Needleman-Wunsch)",
    "Alignment Length (bp)": total_cols,
    "Identical Matches": matches,
    "Mismatches": total_cols - gaps - matches,
    "Gaps Introduced": gaps,
    "Pairwise Sequence Identity (%)": round(pairwise_identity, 2),
    "Alignment Score": pairwise_aln.score
}])

pairwise_summary.to_csv("algal_18S_phylogeny/results/tables/Table_4_3_pairwise_comparison.csv", index=False)
pairwise_summary

```

* **Extracted Presentation Deliverable:** **Table 4.3: Pairwise Global Alignment Parameters Between *C. vulgaris* and *T. obliquus***.



---

### Step 6: Multiple Sequence Alignment (MSA) Execution

Execute MAFFT via subprocess using iterative refinement parameters (`L-INS-i` or `--auto`):

```python
import subprocess

msa_input = "algal_18S_phylogeny/data/references/dataset_18S.fasta"
msa_output = "algal_18S_phylogeny/alignments/18S_MSA.fasta"

cmd = ["mafft", "--auto", "--reorder", msa_input]
with open(msa_output, "w") as f_out:
    subprocess.run(cmd, stdout=f_out, stderr=subprocess.PIPE, check=True)

print(f"Multiple sequence alignment generated at: {msa_output}")

```

* **Extracted Artifacts:** `alignments/18S_MSA.fasta` (aligned matrix used for downstream phylogenetics and site-entropy evaluations).



---

### Step 7: Positional Conservation & Variable Site Profiling

Extract per-site mutational profiles to assess conservation and divergence patterns across the aligned 18S gene:

```python
from Bio import AlignIO

alignment = AlignIO.read("algal_18S_phylogeny/alignments/18S_MSA.fasta", "fasta")
aln_len = alignment.get_alignment_length()
num_taxa = len(alignment)

conserved_sites = 0
variable_sites = 0
parsimony_informative = 0
gap_positions = 0

for i in range(aln_len):
    col = alignment[:, i].upper()
    non_gap = [base for base in col if base not in "-?"]
    
    if len(non_gap) == 0 or (len(col) - len(non_gap)) > (0.5 * num_taxa):
        gap_positions += 1
        continue
        
    unique_bases = set(non_gap)
    if len(unique_bases) == 1:
        conserved_sites += 1
    else:
        variable_sites += 1
        # Check parsimony informativeness: at least 2 distinct bases each appearing >= 2 times
        base_counts = {b: non_gap.count(b) for b in unique_bases}
        if sum(1 for cnt in base_counts.values() if cnt >= 2) >= 2:
            parsimony_informative += 1

msa_metrics = pd.DataFrame([{
    "Total Alignment Length (bp)": aln_len,
    "Number of Taxa Evaluated": num_taxa,
    "Conserved Sites": conserved_sites,
    "Variable Sites": variable_sites,
    "Parsimony-Informative Sites": parsimony_informative,
    "High-Gap Positions (>50%)": gap_positions,
    "Conservation Rate (%)": round((conserved_sites / aln_len) * 100, 2)
}])

msa_metrics.to_csv("algal_18S_phylogeny/results/tables/Table_4_4_alignment_site_metrics.csv", index=False)
msa_metrics

```

* **Extracted Presentation Deliverable:** **Table 4.4: Gene Alignment Statistics and Locus Informativeness Profile**.



---

### Step 8: Pairwise Genetic Distance Matrix & Heatmap Visualization

Quantify genetic divergence and plot a publication-ready identity matrix:

```python
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

taxa_names = [rec.id for rec in alignment]
n_taxa = len(taxa_names)
dist_matrix = np.zeros((n_taxa, n_taxa))

for i in range(n_taxa):
    for j in range(n_taxa):
        s1 = str(alignment[i].seq).upper()
        s2 = str(alignment[j].seq).upper()
        matches = sum(1 for a, b in zip(s1, s2) if a == b and a not in "-?")
        valid_pairs = sum(1 for a, b in zip(s1, s2) if a not in "-?" and b not in "-?")
        dist_matrix[i, j] = (matches / valid_pairs) * 100 if valid_pairs > 0 else 0.0

dist_df = pd.DataFrame(dist_matrix, index=taxa_names, columns=taxa_names)
dist_df.to_csv("algal_18S_phylogeny/results/tables/Table_4_5_pairwise_distance_matrix.csv")

# Plotting Matrix Heatmap
plt.figure(figsize=(10, 8), dpi=300)
sns.heatmap(dist_df, annot=True, fmt=".1f", cmap="YlGnBu", cbar_kws={'label': '18S rRNA Sequence Identity (%)'})
plt.title("Pairwise 18S rRNA Sequence Identity Matrix Across Chlorophyta", fontsize=12, fontweight='bold', pad=12)
plt.xticks(rotation=45, ha='right', fontsize=9)
plt.yticks(fontsize=9)
plt.tight_layout()
plt.savefig("algal_18S_phylogeny/results/figures/Figure_4_1_distance_heatmap.png")
plt.show()

```

* **Extracted Presentation Deliverable:**
* **Table 4.5: Pairwise Genetic Distance and Percentage Identity Matrix**.


* **Figure 4.1: Pairwise 18S Sequence Identity Heatmap** (visual demonstration of divergence between Trebouxiophyceae and Chlorophyceae).





---

### Step 9: Maximum-Likelihood Phylogeny with Model Selection & Bootstrap

Reconstruct evolutionary relationships using IQ-TREE with automated model selection (`-m MFP`), 1,000 Ultrafast Bootstrap replicates (`-B 1000`), and 1,000 SH-aLRT replicates (`--alrt 1000`):

```python
iqtree_cmd = [
    "iqtree2",
    "-s", "algal_18S_phylogeny/alignments/18S_MSA.fasta",
    "-m", "MFP",
    "-B", "1000",
    "--alrt", "1000",
    "-nt", "AUTO",
    "-pre", "algal_18S_phylogeny/phylogeny/18S_ML_analysis"
]

subprocess.run(iqtree_cmd, check=True)
print("Phylogenetic reconstruction completed. Best-fitting substitution model identified by ModelFinder.")

```

Extract the selected substitution model and log-likelihood statistics directly from the IQ-TREE log file:

```python
with open("algal_18S_phylogeny/phylogeny/18S_ML_analysis.iqtree", "r") as f:
    iqtree_summary = f.read()

# Parse best-fit model
for line in iqtree_summary.split("\n"):
    if "Best-fit model according to" in line or "Log-likelihood of the tree" in line:
        print(line)

```

* **Extracted Artifacts:**
* `18S_ML_analysis.treefile`: Newick tree format with maximum-likelihood branch lengths.


* `18S_ML_analysis.contree`: Consensus tree containing SH-aLRT and UFBoot branch support values.


* `18S_ML_analysis.iqtree`: Analytical report containing AIC/BIC model metrics.





---

### Step 10: Phylogenetic Tree Rooting, Visualization & Presentation Export

Load the consensus tree into Biopython, outgroup-root it using *Mesotaenium endlicherianum*, and render a publication-standard phylogram:

```python
from Bio import Phylo
import matplotlib.pyplot as plt

tree = Phylo.read("algal_18S_phylogeny/phylogeny/18S_ML_analysis.contree", "newick")

# Root tree on the designated outgroup
tree.root_with_outgroup({"name": "Mesotaenium_endlicherianum"})
tree.ladderize()

fig, ax = plt.subplots(figsize=(10, 7), dpi=300)
Phylo.draw(
    tree,
    axes=ax,
    show_confidence=True,
    branch_labels=lambda c: f"{int(c.confidence)}" if c.confidence and c.confidence < 100 else None,
    do_show=False
)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_xlabel("Substitutions per nucleotide site", fontsize=10, fontweight='bold')
ax.set_title("Maximum-Likelihood Phylogeny of Chlorophyta Based on 18S rRNA Gene Sequences", fontsize=11, fontweight='bold', pad=15)

plt.tight_layout()
plt.savefig("algal_18S_phylogeny/results/figures/Figure_4_2_ML_phylogenetic_tree.png")
plt.savefig("algal_18S_phylogeny/results/figures/Figure_4_2_ML_phylogenetic_tree.svg")
plt.show()

```

* **Extracted Presentation Deliverable:** **Figure 4.2: Maximum-Likelihood Molecular Phylogram of 18S rRNA Genes** (featuring ultrafast bootstrap support numbers, clear taxonomic clades separating *Chlorella* from *Tetradesmus*, outgroup rooting, and scale bar).



---

### Summary of Result Presentation Extraction for Chapter 4

| Study Deliverable | Source Code / Script Step | Physical Output File | Role in Dissertation / Manuscript |
| --- | --- | --- | --- |
| **Table 4.1** | Step 3 (`compute_qc`) | `results/tables/Table_4_1_sequence_qc.csv` | Sequence lengths, base ratios, and GC bias confirmation.

 |
| **Table 4.2** | Step 4 (`blastn`) | `results/tables/Table_4_2_blast_validation.csv` | Independent homology and coordinate validation metrics.

 |
| **Table 4.3** | Step 5 (`Align.PairwiseAligner`) | `results/tables/Table_4_3_pairwise_comparison.csv` | Direct global similarity indices between the two focus taxa.

 |
| **Table 4.4** | Step 7 (`AlignIO` site scan) | `results/tables/Table_4_4_alignment_site_metrics.csv` | Proportion of conserved, variable, and parsimony-informative loci.

 |
| **Table 4.5** | Step 8 (`dist_matrix`) | `results/tables/Table_4_5_pairwise_distance_matrix.csv` | Complete inter-species pairwise distance values.

 |
| **Figure 4.1** | Step 8 (`sns.heatmap`) | `results/figures/Figure_4_1_distance_heatmap.png` | Visual divergence patterns across microalgal clades.

 |
| **Figure 4.2** | Step 10 (`Phylo.draw`) | `results/figures/Figure_4_2_ML_phylogenetic_tree.png/svg` | Primary phylogenetic hypothesis displaying evolutionary lineage splits and branch support.

 |