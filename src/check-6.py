import pandas as pd

df = pd.read_csv('../data/peptide_conditional_norm.csv')

# Check for sequences with leading/trailing whitespace
mask = df['sequence'] != df['sequence'].str.strip()
bad = df[mask]
print(f"Rows with whitespace in sequence: {len(bad)}")
if len(bad) > 0:
    print(bad[['sequence', 'Activity']].head(20))
    print("\nActual sequence values (repr):")
    for val in bad['sequence'].head(5):
        print(repr(val))
