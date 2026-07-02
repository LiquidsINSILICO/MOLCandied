import pandas as pd
import re

df = pd.read_csv('../data/peptide_conditional_norm.csv')

print(f"Rows before cleaning: {len(df)}")

# Strip all whitespace (spaces, tabs, newlines) from sequence and Activity
df['sequence'] = df['sequence'].str.replace(r'\s+', '', regex=True)
df['Activity'] = df['Activity'].str.replace(r'\s+', ' ', regex=True).str.strip()

# Verify
bad = df[df['sequence'] != df['sequence'].str.strip()]
print(f"Rows with whitespace remaining in sequence: {len(bad)}")

# Also remove any rows where sequence became empty after stripping
before = len(df)
df = df[df['sequence'].str.len() > 0].reset_index(drop=True)
print(f"Rows removed as empty after stripping: {before - len(df)}")

# Validate only standard amino acids remain
AA_SET = set("ACDEFGHIKLMNPQRSTVWY")
invalid = df[~df['sequence'].apply(lambda s: all(c in AA_SET for c in s))]
print(f"Rows with non-standard amino acids: {len(invalid)}")
if len(invalid) > 0:
    print(invalid[['sequence', 'Activity']].head(10))

# Save cleaned file
df.to_csv('../data/peptide_conditional_norm-2.csv', index=False)
print(f"\nSaved cleaned file with {len(df)} rows")

# Spot check
print("\nSample sequences after cleaning:")
print(df['sequence'].head(10).tolist())
print("\nSample Activity after cleaning:")
print(df['Activity'].head(10).tolist())
