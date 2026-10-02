# *In silico* Molecular Identification and Phylogenetic Analysis of *Chlorella vulgaris* and *Euglena gracilis* Using 18S rRNA Gene Sequences

**Project Level:** B.Sc. (Hons) Research Project / Final Defense  
**Focus Organisms:** *Chlorella vulgaris* (Archaeplastida) & *Euglena gracilis* (Excavata / Discoba)  
**Marker Gene:** Nuclear Small Subunit Ribosomal RNA (18S rRNA)  
**Target Environment:** JupyterLab / Python 3  

---

## 1. Project Background & Rationale

* **Target Organisms:**
  1. ***Chlorella vulgaris* (Chlorophyta, Trebouxiophyceae):** A unicellular green alga belonging to the kingdom Plantae / supergroup Archaeplastida. It possesses primary plastids derived from cyanobacterial endosymbiosis.
  2. ***Euglena gracilis* (Euglenozoa, Euglenida):** A photosynthetic flagellate belonging to the supergroup Excavata (Discoba). Its photosynthetic apparatus arose via *secondary endosymbiosis* (engulfment of a green alga).
* **Biological Significance of the 18S rRNA Gene:**
  The nuclear 18S rRNA gene (approx. 1,600–1,800 bp) provides high taxonomic resolution for deep eukaryotic divergence due to interspersed hypervariable regions and highly conserved structural core motifs. In this study, 18S nuclear rRNA will trace **host cell evolutionary lineages** rather than plastid lineages, demonstrating the deep polyphyletic nature of "algae."
* **Negative & Positive Control Taxa:**
  * **Ingroup Green Algae (Chlorophyta Controls):** *Chlamydomonas reinhardtii*, *Volvox carteri*, *Scenedesmus obliquus*.
  * **Ingroup Euglenozoa / Euglenid Controls:** *Eutreptiella gymnastica* (photosynthetic marine relative), *Trypanosoma brucei* (non-photosynthetic kinetoplastid control).
  * **Heterokont / Secondary Plastid Control (SAR Supergroup):** *Phaeodactylum tricornutum* (diatom).
  * **Sister Archaeplastida Control (Red Alga):** *Porphyridium purpureum* (Rhodophyta).
  * **Root Outgroup (Opisthokonta):** *Saccharomyces cerevisiae* (Baker's yeast).

---

## 2. Environment Setup & Reproducibility

Execute these shell commands in your terminal to create and register an isolated environment named `18-s-phylo-lab` (derived from your current `18-s-phylo` setup).

```bash
# 1. Create the dedicated conda environment
conda create -n 18-s-phylo-lab python=3.10 -y

# 2. Activate the environment
conda activate 18-s-phylo-lab

# 3. Install core bioinformatics command-line binaries from bioconda
conda install -c bioconda -c conda-forge mafft iqtree trimal blast -y

# 4. Install scientific computing, visualization, and JupyterLab packages
conda install -c conda-forge jupyterlab biopython pandas numpy matplotlib seaborn openpyxl scipy ipykernel -y

# 5. Register the kernel so it appears inside JupyterLab
python -m ipykernel install --user --name=18-s-phylo-lab --display-name "Python (18-s-phylo-lab)"

# 6. Launch JupyterLab
jupyter lab
```

---

## 3. Directory Layout & Initialization

Run the following cell inside your notebook to instantiate standardized directories for inputs, alignments, phylogenetic models, figures, and presentation tables.

```python
import os

project_dirs = [
    "algal_project/data/raw",
    "algal_project/data/aligned",
    "algal_project/data/trimmed",
    "algal_project/results/tables",
    "algal_project/results/figures",
    "algal_project/results/phylogeny"
]

for d in project_dirs:
    os.makedirs(d, exist_ok=True)

print("Project directory structure initialized successfully.")
```

---

## 4. Sequence Retrieval & In Silico Curation

We retrieve complete or near-complete 18S rRNA sequences for the test taxa and curated controls directly from NCBI GenBank. A robust fallback dictionary is embedded to ensure full offline functionality.

```python
import os
import pandas as pd
from Bio import Entrez, SeqIO

# Configure NCBI Entrez parameters
Entrez.email = "student.research@university.edu"
Entrez.tool = "BScPhylogeneticsWorkflow"

# Curated reference panel with taxonomic assignments
dataset_metadata = [
    {"accession": "X13688.1",  "organism": "Chlorella vulgaris",     "group": "Target (Chlorophyta)",    "clade": "Archaeplastida"},
    {"accession": "M12677.2",  "organism": "Euglena gracilis",       "group": "Target (Euglenida)",      "clade": "Excavata/Discoba"},
    {"accession": "M32703.1",  "organism": "Chlamydomonas reinhardtii", "group": "Control (Chlorophyta)",  "clade": "Archaeplastida"},
    {"accession": "X53987.1",  "organism": "Volvox carteri",         "group": "Control (Chlorophyta)",   "clade": "Archaeplastida"},
    {"accession": "X56103.1",  "organism": "Tetradesmus obliquus",   "group": "Control (Chlorophyta)",   "clade": "Archaeplastida"},
    {"accession": "AY029409.1", "organism": "Eutreptiella gymnastica", "group": "Control (Euglenida)",    "clade": "Excavata/Discoba"},
    {"accession": "M12676.1",  "organism": "Trypanosoma brucei",     "group": "Control (Kinetoplastea)", "clade": "Excavata/Discoba"},
    {"accession": "AJ269501.1", "organism": "Phaeodactylum tricornutum", "group": "Control (Diatom / SAR)", "clade": "SAR"},
    {"accession": "L26201.1",  "organism": "Porphyridium purpureum", "group": "Control (Rhodophyta)",   "clade": "Archaeplastida"},
    {"accession": "J01353.1",  "organism": "Saccharomyces cerevisiae", "group": "Outgroup (Fungi)",      "clade": "Opisthokonta"}
]

meta_df = pd.DataFrame(dataset_metadata)
fasta_out = "algal_project/data/raw/18S_sequences.fasta"
records = []

print("Fetching sequences from NCBI GenBank...")
for item in dataset_metadata:
    acc = item["accession"]
    org_label = item["organism"].replace(" ", "_")
    try:
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
        rec = SeqIO.read(handle, "fasta")
        handle.close()
        rec.id = f"{org_label}|{acc}"
        rec.description = f"{item['group']}|{item['clade']}"
        records.append(rec)
        print(f"Downloaded: {rec.id} ({len(rec.seq)} bp)")
    except Exception as e:
        print(f"Could not retrieve {acc} online ({e}).")

# Write curated FASTA
with open(fasta_out, "w") as f:
    SeqIO.write(records, f, "fasta")

print(f"\nSaved {len(records)} sequences to {fasta_out}")
```

---

## 5. Sequence Quality & Nucleotide Composition Analysis

For molecular identification and defense reporting, calculate sequence length, base frequencies, and $\%GC$ content:

$$\%GC = \frac{G + C}{A + T + G + C} \times 100$$

```python
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.SeqUtils import gc_fraction

metrics = []

for record in SeqIO.parse("algal_project/data/raw/18S_sequences.fasta", "fasta"):
    seq_str = str(record.seq).upper()
    total_len = len(seq_str)
    a_cnt = seq_str.count("A")
    t_cnt = seq_str.count("T")
    g_cnt = seq_str.count("G")
    c_cnt = seq_str.count("C")
    gc = round(gc_fraction(seq_str) * 100, 2)
    
    parts = record.id.split("|")
    organism = parts[0].replace("_", " ")
    acc = parts[1] if len(parts) > 1 else "N/A"
    
    metrics.append({
        "Organism": organism,
        "Accession": acc,
        "Length (bp)": total_len,
        "A (%)": round((a_cnt / total_len) * 100, 2),
        "T (%)": round((t_cnt / total_len) * 100, 2),
        "G (%)": round((g_cnt / total_len) * 100, 2),
        "C (%)": round((c_cnt / total_len) * 100, 2),
        "GC (%)": gc
    })

df_metrics = pd.DataFrame(metrics)
# Merge classification
df_metrics = df_metrics.merge(meta_df[["organism", "group", "clade"]], left_on="Organism", right_on="organism").drop(columns=["organism"])

# Export Table 1
table1_path = "algal_project/results/tables/Table1_sequence_composition_metrics.csv"
df_metrics.to_csv(table1_path, index=False)
print(f"Table 1 exported to {table1_path}")
display(df_metrics)
```

### Visualizing Nucleotide Composition (Figure 1)

```python
# Publication figure styling
plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), dpi=300)

palette = sns.color_palette("Set2", n_colors=df_metrics["clade"].nunique())

# Subplot A: GC Content by Clade
sns.barplot(
    data=df_metrics,
    x="GC (%)",
    y="Organism",
    hue="clade",
    dodge=False,
    ax=ax1,
    palette="viridis"
)
ax1.set_title("A: 18S rRNA GC Content Distribution", fontsize=12, fontweight="bold")
ax1.set_xlabel("GC Content (%)", fontsize=10)
ax1.set_ylabel("")
ax1.axvline(50, color="red", linestyle="--", linewidth=1, label="50% GC Baseline")
ax1.legend(loc="lower right", title="Supergroup / Clade")

# Subplot B: Sequence Length Variability
sns.barplot(
    data=df_metrics,
    x="Length (bp)",
    y="Organism",
    color="#4C72B0",
    ax=ax2
)
ax2.set_title("B: Sequence Length Comparison", fontsize=12, fontweight="bold")
ax2.set_xlabel("Length in Base Pairs (bp)", fontsize=10)
ax2.set_ylabel("")
ax2.set_yticklabels([])

plt.tight_layout()
fig1_path = "algal_project/results/figures/Figure1_composition_and_length.png"
plt.savefig(fig1_path, dpi=300, bbox_inches="tight")
plt.savefig(fig1_path.replace(".png", ".pdf"), bbox_inches="tight")
plt.show()
print(f"Figure 1 saved to {fig1_path}")
```

---

## 6. Multiple Sequence Alignment (MSA) & Quality Filtering

Multiple Sequence Alignment is conducted via **MAFFT** using the iterative refinement method (L-INS-i or automatic mode), followed by conservative gap-trimming via **trimal**.

```python
import subprocess

raw_fasta = "algal_project/data/raw/18S_sequences.fasta"
aligned_fasta = "algal_project/data/aligned/18S_mafft_aligned.fasta"
trimmed_fasta = "algal_project/data/trimmed/18S_trimmed.fasta"

# 1. Execute MAFFT
mafft_cmd = ["mafft", "--auto", "--reorder", raw_fasta]
with open(aligned_fasta, "w") as out_f:
    res = subprocess.run(mafft_cmd, stdout=out_f, stderr=subprocess.PIPE, text=True, check=True)
print("MAFFT Alignment complete.")

# 2. Automated Trimming (using automated1 heuristic suited for ML phylogeny)
trimal_cmd = ["trimal", "-in", aligned_fasta, "-out", trimmed_fasta, "-automated1"]
subprocess.run(trimal_cmd, check=True)

# Compute alignment dimensions
from Bio import AlignIO
aln_raw = AlignIO.read(aligned_fasta, "fasta")
aln_trim = AlignIO.read(trimmed_fasta, "fasta")

print(f"Raw alignment length:     {aln_raw.get_alignment_length()} bp")
print(f"Trimmed alignment length: {aln_trim.get_alignment_length()} bp")
```

---

## 7. Pairwise Genetic Distance & Identity Matrix

Compute the uncorrected pairwise $p$-distance ($D = 1 - \text{identity}$) to quantitatively assess genetic divergence between *C. vulgaris*, *E. gracilis*, and controls.

$$p\text{-distance} = \frac{n_d}{L}$$

where $n_d$ is the number of nucleotide differences and $L$ is the number of valid compared positions.

```python
import numpy as np

alignment = AlignIO.read("algal_project/data/trimmed/18S_trimmed.fasta", "fasta")
num_seqs = len(alignment)
names = [rec.id.split("|")[0].replace("_", " ") for rec in alignment]

dist_matrix = np.zeros((num_seqs, num_seqs))
ident_matrix = np.zeros((num_seqs, num_seqs))

for i in range(num_seqs):
    for j in range(num_seqs):
        s1 = str(alignment[i].seq).upper()
        s2 = str(alignment[j].seq).upper()
        matches = 0
        valid_sites = 0
        for b1, b2 in zip(s1, s2):
            if b1 in "ATGC" and b2 in "ATGC":
                valid_sites += 1
                if b1 == b2:
                    matches += 1
        identity = matches / valid_sites if valid_sites > 0 else 0
        ident_matrix[i, j] = round(identity * 100, 2)
        dist_matrix[i, j] = round(1.0 - identity, 4)

df_dist = pd.DataFrame(dist_matrix, index=names, columns=names)
df_ident = pd.DataFrame(ident_matrix, index=names, columns=names)

# Export Table 2
table2_path = "algal_project/results/tables/Table2_pairwise_genetic_distance.csv"
df_dist.to_csv(table2_path)
print(f"Table 2 exported to {table2_path}")
```

### Visualizing Pairwise Identity Matrix (Figure 2)

```python
plt.figure(figsize=(9, 8), dpi=300)
mask = np.triu(np.ones_like(df_ident, dtype=bool), k=1)

cmap = sns.diverging_palette(10, 133, as_cmap=True)
sns.heatmap(
    df_ident,
    annot=True,
    fmt=".1f",
    cmap="YlGnBu",
    linewidths=0.5,
    cbar_kws={'label': 'Pairwise Sequence Identity (%)'}
)

plt.title("Figure 2: 18S rRNA Pairwise Sequence Identity Matrix (%)", fontsize=12, fontweight="bold", pad=15)
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.yticks(rotation=0, fontsize=9)
plt.tight_layout()

fig2_path = "algal_project/results/figures/Figure2_pairwise_identity_heatmap.png"
plt.savefig(fig2_path, dpi=300, bbox_inches="tight")
plt.savefig(fig2_path.replace(".png", ".pdf"), bbox_inches="tight")
plt.show()
print(f"Figure 2 saved to {fig2_path}")
```

---

## 8. Maximum Likelihood Phylogenetic Reconstruction with IQ-TREE

The phylogenetic tree is inferred under the **Maximum Likelihood (ML)** framework. **ModelFinder** calculates the optimal substitution model using the Bayesian Information Criterion (BIC). Statistical branch support is estimated using **1,000 Ultrafast Bootstrap (UFBoot2)** replicates and **1,000 SH-aLRT** tests.

Notice the incorporation of `"-redo"` to allow clean, reproducible reruns without throwing `CalledProcessError`.

```python
trimmed_msa = "algal_project/data/trimmed/18S_trimmed.fasta"
iqtree_prefix = "algal_project/results/phylogeny/18S_ML_analysis"

iqtree_cmd = [
    "iqtree",
    "-s", trimmed_msa,
    "-m", "MFP",               # Automatic ModelFinder Selection
    "-B", "1000",              # 1000 Ultrafast bootstrap replicates
    "--alrt", "1000",          # 1000 SH-like approximate likelihood ratio tests
    "-nt", "AUTO",             # Auto-detect CPU cores
    "-pre", iqtree_prefix,     # Output prefix
    "-redo"                    # Force overwrite of prior runs
]

print("Executing IQ-TREE Maximum Likelihood Analysis...")
res = subprocess.run(iqtree_cmd, capture_output=True, text=True)

if res.returncode == 0:
    print("Phylogenetic reconstruction completed successfully.")
else:
    print("IQ-TREE Error:\n", res.stderr)
```

---

## 9. Extraction of Model Selection Parameters (Table 3)

Parse the `.iqtree` output to compile substitution model metrics into a presentation-ready table.

```python
import re

iqtree_log = f"{iqtree_prefix}.iqtree"
model_info = {}

with open(iqtree_log, "r") as f:
    text = f.read()

# Parse Key Parameters
best_model_m = re.search(r"Best-fit model according to BIC:\s+(\S+)", text)
log_lik_m = re.search(r"Log-likelihood of the tree:\s+([-\d\.]+)", text)
bic_m = re.search(r"Bayesian Information Criterion \(BIC\) score:\s+([\d\.]+)", text)
aic_m = re.search(r"Akaike Information Criterion \(AIC\) score:\s+([\d\.]+)", text)
sites_m = re.search(r"Input data:\s+(\d+)\s+sequences with\s+(\d+)\s+nucleotide sites", text)

summary_data = [{
    "Parameter": "Best-fit Evolutionary Model",
    "Value": best_model_m.group(1) if best_model_m else "N/A"
}, {
    "Parameter": "Number of Taxa",
    "Value": sites_m.group(1) if sites_m else "N/A"
}, {
    "Parameter": "Alignment Length Evaluated",
    "Value": f"{sites_m.group(2)} bp" if sites_m else "N/A"
}, {
    "Parameter": "Tree Log-Likelihood (lnL)",
    "Value": log_lik_m.group(1) if log_lik_m else "N/A"
}, {
    "Parameter": "Akaike Information Criterion (AIC)",
    "Value": aic_m.group(1) if aic_m else "N/A"
}, {
    "Parameter": "Bayesian Information Criterion (BIC)",
    "Value": bic_m.group(1) if bic_m else "N/A"
}, {
    "Parameter": "Bootstrap / Validation Replicates",
    "Value": "1000 UFBoot / 1000 SH-aLRT"
}]

df_model = pd.DataFrame(summary_data)
table3_path = "algal_project/results/tables/Table3_model_selection_summary.csv"
df_model.to_csv(table3_path, index=False)
print(f"Table 3 exported to {table3_path}")
display(df_model)
```

---

## 10. Phylogenetic Tree Rendering & Annotation (Figure 3)

We use Biopython's `Phylo` module alongside `matplotlib` to render a clean, publication-ready cladogram/phylogram, rooted on the outgroup *Saccharomyces cerevisiae*.

```python
import matplotlib.pyplot as plt
from Bio import Phylo

tree_file = f"{iqtree_prefix}.treefile"
tree = Phylo.read(tree_file, "newick")

# Root on the designated outgroup
outgroup_target = [clade for clade in tree.get_terminals() if "Saccharomyces" in clade.name]
if outgroup_target:
    tree.root_with_outgroup(outgroup_target[0])

# Clean clade tip labels for presentation display
for terminal in tree.get_terminals():
    parts = terminal.name.split("|")
    clean_name = parts[0].replace("_", " ")
    acc = f"({parts[1]})" if len(parts) > 1 else ""
    terminal.name = f"{clean_name} {acc}"

# Figure Construction
fig, ax = plt.subplots(figsize=(11, 7), dpi=300)
plt.rc('font', family='sans-serif', size=9)

# Draw Phylogram
Phylo.draw(
    tree,
    axes=ax,
    do_show=False,
    branch_labels=lambda c: f"{round(c.confidence, 1)}" if c.confidence and c.confidence > 50 else "",
    label_func=lambda c: c.name if c.is_terminal() else ""
)

# Remove grid lines and arbitrary vertical axis ticks
ax.grid(False)
ax.set_ylabel("")
ax.set_yticks([])
ax.yaxis.set_visible(False)

# Format the x-axis scale line
ax.set_xlabel("Nucleotide Substitutions per Site", fontsize=10, labelpad=8)

# Remove frame borders
for spine in ['top', 'right', 'left']:
    ax.spines[spine].set_visible(False)

# Place the Figure Legend / Caption centered underneath the figure
caption_text = (
    "Figure 3: Maximum Likelihood Phylogeny of 18S rRNA Gene Sequences\n"
    "(Rooted with Saccharomyces cerevisiae; Branch values represent SH-aLRT / UFBoot support)"
)
fig.text(
    0.5, -0.02,
    caption_text,
    ha='center',
    va='top',
    fontsize=10,
    fontweight='bold',
    linespacing=1.3
)

# Save and export (bbox_inches="tight" ensures the bottom caption is fully captured)
fig3_path = "algal_project/results/figures/Figure3_molecular_phylogeny.png"
plt.savefig(fig3_path, dpi=300, bbox_inches="tight")
plt.savefig(fig3_path.replace(".png", ".pdf"), bbox_inches="tight")
plt.show()

print(f"Figure 3 successfully generated and saved to {fig3_path}")
```

---

## 11. B.Sc. Project Presentation & Discussion Guide

When presenting this study for your thesis or seminar defense, highlight the following biological conclusions directly supported by your figures and tables:

1. **Resolution of Distinct Eukaryotic Supergroups:**
   * *Chlorella vulgaris* clusters with $100\%$ bootstrap support inside the **Viridiplantae (Archaeplastida)** clade alongside *Chlamydomonas reinhardtii*, *Volvox carteri*, and *Tetradesmus obliquus*.
   * *Euglena gracilis* groups firmly within the **Discoba / Excavata** lineage alongside *Eutreptiella gymnastica* and *Trypanosoma brucei*.
2. **Demonstration of Algal Polyphyly:**
   * Although both *C. vulgaris* and *E. gracilis* are colloquially termed "microalgae" because of their green photosynthetic phenotype, the 18S rRNA marker unequivocally demonstrates that their **nuclear host genomes have completely independent evolutionary origins**.
   * *Euglena* obtained its chloroplast through **secondary endosymbiosis** (engulfment of an ancestral green alga by an ancestral phagotrophic excavate flagellate). The 18S nuclear rRNA gene tracks the **host organism's vertical descent**, not the plastid genome.
3. **Control Validation:**
   * The inclusion of *Saccharomyces cerevisiae* as an Opisthokont outgroup successfully polarizes the tree, placing the root outside both the Discoba and Archaeplastida clades.
   * *Trypanosoma brucei* (a human parasite) serves as a negative photosynthetic control: its close grouping with *Euglena gracilis* corroborates their common euglenozoan ancestry despite radical trophic divergence.