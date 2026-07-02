import pandas as pd
import numpy as np
import pickle

from scipy.stats import skew
#from rdkit import Chem
#from rdkit.Chem import Descriptors
from sklearn.preprocessing import StandardScaler
from Bio.SeqUtils import molecular_weight
from Bio.Seq import Seq

# =====================================================
# INPUT / OUTPUT
# =====================================================

INPUT_FILE = "peptides_filtered-conditional.csv"

UNORM_FILE = "peptide_conditional_unorm.csv"
NORM_FILE = "peptide_conditional_norm.csv"

SCALER_FILE = "property_scaler.pkl"


# =====================================================
# KYTE-DOOLITTLE SCALE
# =====================================================

KD = {
    'A': 1.8,
    'R': -4.5,
    'N': -3.5,
    'D': -3.5,
    'C': 2.5,
    'Q': -3.5,
    'E': -3.5,
    'G': -0.4,
    'H': -3.2,
    'I': 4.5,
    'L': 3.8,
    'K': -3.9,
    'M': 1.9,
    'F': 2.8,
    'P': -1.6,
    'S': -0.8,
    'T': -0.7,
    'W': -0.9,
    'Y': -1.3,
    'V': 4.2
}

HYDROPHOBIC = {
    "A",
    "V",
    "I",
    "L",
    "M",
    "F",
    "W",
    "Y",
    "C"
    }

# =====================================================
# SMILES FROM SEQUENCE
# =====================================================
# Exact accepted activity strings for each taste
EXACT_TASTE_MAP = {
    "sweet":  "sweet",
    "bitter": "bitter",
    "umami":  "umami",
    "sour":   "sour",
    "salty":  "salty"
}

def encode_taste(activity):
    activity = str(activity).lower().strip()
    
    # Split by comma or slash to handle multi-taste entries like "sweet, umami"
    parts = {chunk.strip() for chunk in activity.replace("/", ",").split(",")}
    
    return pd.Series({
        taste: int(exact_term in parts)
        for taste, exact_term in EXACT_TASTE_MAP.items()
    })

# =====================================================
# PEPTIDE LENGTH
# =====================================================

def peptide_length(seq):

    return len(str(seq).strip())

# =====================================================
# NET CHARGE
# =====================================================

def peptide_charge(seq):

    seq = str(seq).upper()

    charge = (
          seq.count("K")
        + seq.count("R")
        + 0.1 * seq.count("H")
        - seq.count("D")
        - seq.count("E")
    )

    return charge


# =====================================================
# AROMATIC FRACTION
# =====================================================

def aromatic_fraction(seq):

    seq = str(seq).upper()

    aromatic = (
          seq.count("F")
        + seq.count("W")
        + seq.count("Y")
    )

    if len(seq) == 0:
        return 0

    return aromatic / len(seq)


# =====================================================
# GRAVY
# =====================================================

def gravy(seq):

    seq = str(seq).upper()

    values = [
        KD[x]
        for x in seq
        if x in KD
    ]

    if len(values) == 0:
        return 0

    return np.mean(values)

# =====================================================
# Fraction of HYDROPHOBIC RESIDUES
# =====================================================

def hydrophobic_fraction(seq):

    seq = str(seq).upper()

    if len(seq) == 0:
        return 0

    count = sum(
            aa in HYDROPHOBIC
            for aa in seq
            )

    return count / len(seq)

# =====================================================
# MW FROM Peptides
# =====================================================

from Bio.SeqUtils.ProtParam import ProteinAnalysis

def peptide_mw(seq):
    try:
        seq_clean = str(seq).strip().upper()
        if len(seq_clean) == 0:
            return np.nan
        return ProteinAnalysis(seq_clean).molecular_weight()
    except Exception:
        return np.nan

# =====================================================
# LOAD DATA
# =====================================================

print("="*60)
print("Loading dataset")
print("="*60)

df = pd.read_csv(INPUT_FILE)

print("Rows:", len(df))

print("\nColumns:")
print(df.columns.tolist())

# =====================================================
# CLEAN RAW INPUT IMMEDIATELY
# =====================================================
print("\nCleaning raw sequences...")

# Remove all whitespace (tabs, spaces) from sequence
df['sequence'] = df['sequence'].str.replace(r'\s+', '', regex=True)

# Uppercase all sequences
df['sequence'] = df['sequence'].str.upper()

# Clean Activity column — strip whitespace
df['Activity'] = df['Activity'].str.replace(r'\s+', ' ', regex=True).str.strip()

# Drop empty sequences
before = len(df)
df = df[df['sequence'].str.len() > 0].reset_index(drop=True)
print(f"Dropped {before - len(df)} empty sequences")

# .upper() inside lambda is a safety net only — sequence is already uppercased above
AA_SET = set("ACDEFGHIKLMNPQRSTVWY")
before = len(df)
df = df[df['sequence'].apply(
    lambda s: all(c in AA_SET for c in s.upper())
)].reset_index(drop=True)
print(f"Dropped {before - len(df)} sequences with non-standard amino acids")

print(f"Clean sequences remaining: {len(df)}")

# Spot check
print("\nSample sequences after cleaning:")
print(df['sequence'].head(10).tolist())
print("Sample Activity after cleaning:")
print(df['Activity'].head(10).tolist())
#================================================================================
# DESCRIPTORS
# ===============================================================================

print("\nCalculating descriptors...")

df["length"] = df["sequence"].apply(
    peptide_length
)

df["charge"] = df["sequence"].apply(
    peptide_charge
)

df["aromatic"] = df["sequence"].apply(
    aromatic_fraction
)

df["gravy"] = df["sequence"].apply(
    gravy
)

df["hydrophobic_fraction"] = df["sequence"].apply(
        hydrophobic_fraction
        )

#df["mw"] = df["smiles"].apply(
#    molecular_weight
#)
df["mw"] = df["sequence"].apply(
        peptide_mw
        )

before = len(df)
df = df.dropna(subset=["mw"]).reset_index(drop=True)
print(f"Dropped {before - len(df)} rows where MW calculation failed")

# =====================================================
# TASTE LABELS
# =====================================================

print("Encoding taste labels...")

taste_df = df["Activity"].apply(
    encode_taste
)

df = pd.concat(
    [df, taste_df],
    axis=1
)

# =====================================================
# DATA SPLIT
# =====================================================
df = df.sample(
        frac=1,
        random_state=42
        ).reset_index(drop=True)

n = len(df)

train_end = int(0.80*n)
val_end   = int(0.90*n)

df["source"] = "train"

df.loc[
    train_end:val_end-1,
    "source"
] = "val"

df.loc[
    val_end:,
    "source"
] = "test"

print("\nDataset split")
print("="*60)

print(df["source"].value_counts())
# =====================================================
# SCAFFOLD
# =====================================================

#df["scaffold_smiles"] = df["smiles"]

# =====================================================
# SAVE UNNORMALIZED
# =====================================================

unorm_cols = [

    "sequence",
    "Activity",
#    "smiles",
    "source",
#    "scaffold_smiles",

    "sweet",
    "bitter",
    "umami",
    "sour",
    "salty",

    "length",
    "charge",
    "mw",
    "aromatic",
    "gravy",
    "hydrophobic_fraction"
]

df_unorm = df[unorm_cols]

df_unorm.to_csv(
    UNORM_FILE,
    index=False
)

print("\nSaved:", UNORM_FILE)


# =====================================================
# RAW STATISTICS
# =====================================================

print("\nRaw statistics")
print("="*60)

for col in [
    "length",
    "charge",
    "mw",
    "aromatic",
    "gravy",
    "hydrophobic_fraction"
]:

    print(
        f"{col:12s}"
        f" mean={df[col].mean():8.3f}"
        f" std={df[col].std():8.3f}"
        f" min={df[col].min():8.3f}"
        f" max={df[col].max():8.3f}"
    )


# =====================================================
# NORMALIZATION
# =====================================================

continuous = [
    "length",
    "charge",
    "mw",
    "aromatic",
    "gravy",
    "hydrophobic_fraction"
]

print("\nskewness analysis")
print("="*60)

transform_cols = []
skew_report = []

for col in continuous:
    sk = skew(df[col].dropna())
    skew_report.append({
        "feature": col,
        "skewness": sk
        })

    print(f"{col:12s} skew = {sk:8.3f}")
    if sk > 1 and col in ["length", "mw"]:
        transform_cols.append(col)

# -----------------------------------------------------
# SAVE SKEWNESS REPORT
# -----------------------------------------------------

pd.DataFrame(skew_report).to_csv(
    "descriptor_skewness.csv",
    index=False
)

print("\nSaved: descriptor_skewness.csv")

#==================================================
# Log transform
#==================================================

print("\nApplying log transform where needed")
print("="*60)
for col in transform_cols:
    print(f"Log-transforming {col}")
    df[f"{col}_raw"] = df[col]
    df[col] = np.log(df[col])

print("\nStatistics after log transform")
print("="*60)

for col in continuous:

    print(
            f"{col:12s}"
            f" mean={df[col].mean():8.3f}"
            f" std={df[col].std():8.3f}"
            f" min={df[col].min():8.3f}"
            f" max={df[col].max():8.3f}"
            ) 

print("\nApplying Z-score normalization...")
print("="*60)

scaler = StandardScaler()

df_norm = df.copy()

train_mask = df_norm["source"] == "train"

scaler.fit(
        df_norm.loc[train_mask, continuous]
        )

df_norm[continuous] = scaler.transform(
    df_norm[continuous]
)

#==========================================================
import matplotlib.pyplot as plt

print("\nGenerating descriptor plots")

for feat in continuous:

    fig, ax = plt.subplots(
        1,
        2,
        figsize=(10,4)
    )

    raw_col = feat+"_raw" if feat+"_raw" in df.columns else feat

    ax[0].hist(
            df[raw_col],
            bins=30,
            density=True
            )

    ax[0].set_title(
        f"Raw {feat}"
    )

    ax[1].hist(
        df_norm[feat],
        bins=30,
        density=True
    )

    ax[1].set_title(
        f"Normalized {feat}"
    )

    plt.tight_layout()

    plt.savefig(
        f"{feat}_raw_vs_norm.png",
        dpi=300
    )

    plt.close()
#==========================================================
taste_cols = [
    "sweet",
    "bitter",
    "umami",
    "sour",
    "salty"
]

counts = [
    df[t].sum()
    for t in taste_cols
]

plt.figure(figsize=(6,4))

plt.bar(
    taste_cols,
    counts
)

for i,v in enumerate(counts):
    plt.text(i, v, str(int(v)), ha="center")

plt.ylabel(
    "Number of peptides"
)

plt.title(
    "Taste Distribution"
)

plt.tight_layout()

plt.savefig(
    "taste_distribution.png",
    dpi=300
)

plt.close()
#==========================================================
from collections import Counter

AA = list(
    "ACDEFGHIKLMNPQRSTVWY"
)

fig, axes = plt.subplots(
    5,
    1,
    figsize=(12,18)
)

for i,taste in enumerate(taste_cols):

    subset = df[
        df[taste] == 1
    ]

    counter = Counter()

    for seq in subset["sequence"]:
        counter.update(
            str(seq)
        )

    total = sum(
        counter.values()
    )

    percentages = [
        100*counter[a]/total
        if total > 0 else 0
        for a in AA
    ]

    axes[i].bar(
        AA,
        percentages
    )

    axes[i].set_title(
        f"{taste.capitalize()} peptides"
    )

    axes[i].set_ylabel(
        "%"
    )

plt.tight_layout()

plt.savefig(
    "AA_composition_all_tastes.png",
    dpi=300
)

plt.close()
#==========================================================
global_counter = Counter()

for seq in df["sequence"]:
    global_counter.update(
        str(seq)
    )

global_total = sum(
    global_counter.values()
)

global_freq = {
    aa:
    global_counter[aa]/global_total
    for aa in AA
}

fig, axes = plt.subplots(
    5,
    1,
    figsize=(12,18)
)

for i,taste in enumerate(taste_cols):

    subset = df[
        df[taste] == 1
    ]

    counter = Counter()

    for seq in subset["sequence"]:
        counter.update(
            str(seq)
        )

    total = sum(
        counter.values()
    )

    enrichment = []

    for aa in AA:

        class_freq = (
            counter[aa]/total
            if total > 0 else 0
        )

        enrich = (
            class_freq /
            global_freq[aa]
            if global_freq[aa] > 0
            else 0
        )

        enrichment.append(
            enrich
        )

    axes[i].bar(
        AA,
        enrichment
    )

    axes[i].axhline(
        1,
        linestyle="--"
    )

    axes[i].set_title(
        f"{taste.capitalize()} enrichment"
    )

plt.tight_layout()

plt.savefig(
    "AA_enrichment_all_tastes.png",
    dpi=300
)

plt.close()
#==========================================================
# =====================================================
# SAVE SCALER
# =====================================================

with open(
    SCALER_FILE,
    "wb"
) as f:

    pickle.dump(
        scaler,
        f
    )

print("\nSaved scaler:", SCALER_FILE)

# =====================================================
# SAVE DESCRIPTOR STATS
# =====================================================

stats = []

for col in continuous:

    stats.append({
        "feature": col,
        "mean": df[col].mean(),
        "std": df[col].std(),
        "min": df[col].min(),
        "max": df[col].max()
    })

pd.DataFrame(stats).to_csv(
    "descriptor_statistics.csv",
    index=False
)

print("\nSaved: descriptor_statistics.csv")

# =====================================================
# SAVE NORMALIZED DATASET
# =====================================================

df_norm[unorm_cols].to_csv(
    NORM_FILE,
    index=False
)

print("Saved:", NORM_FILE)


# =====================================================
# NORMALIZED STATISTICS
# =====================================================

print("\nNormalized statistics")
print("="*60)

for col in continuous:

    print(
        f"{col:12s}"
        f" mean={df_norm[col].mean():8.3f}"
        f" std={df_norm[col].std():8.3f}"
        f" min={df_norm[col].min():8.3f}"
        f" max={df_norm[col].max():8.3f}"
    )


# =====================================================
# PREVIEW
# =====================================================

print("\nFirst 5 rows (unnormalized)")
print(
    df_unorm.head()
)

print("\nFirst 5 rows (normalized)")
print(
    df_norm[unorm_cols].head()
)

print("\nFinal dataset size")
print("="*60)

print(f"Total peptides : {len(df)}")

for t in [
    "sweet",
    "bitter",
    "umami",
    "sour",
    "salty"
]:
    print(
        f"{t:10s} : {int(df[t].sum())}"
    )

print("\nDone.")
