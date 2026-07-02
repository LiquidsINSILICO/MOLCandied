import pandas as pd

df = pd.read_csv('../data/peptide_conditional_norm.csv')

properties = ['sweet', 'bitter', 'umami', 'sour', 'salty', 
              'length', 'charge', 'mw', 'aromatic', 'gravy', 'hydrophobic_fraction']

print("Total rows:", len(df))

for col in properties:
    numeric_check = pd.to_numeric(df[col], errors='coerce')
    bad_rows = df[numeric_check.isna() & df[col].notna()]
    if len(bad_rows) > 0:
        print(f"\n'{col}' has {len(bad_rows)} non-numeric value(s):")
        print(bad_rows[['sequence', 'Activity', col]])
    
    nan_rows = df[col].isna().sum()
    if nan_rows > 0:
        print(f"'{col}' has {nan_rows} NaN (missing) value(s)")
