import pandas as pd

df = pd.read_csv('../data/peptide_conditional_norm.csv')

properties = ['sweet','bitter','umami','sour','salty','length',
              'charge','mw','aromatic','gravy','hydrophobic_fraction']

print("Checking for string 'nan' values hiding in numeric columns...")
found_any = False

for col in properties:
    for idx, val in enumerate(df[col]):
        # Check if the actual Python type is string, not float/int
        if isinstance(val, str):
            print(f"FOUND IT: column='{col}' row={idx} value='{val}' sequence='{df['sequence'].iloc[idx]}'")
            found_any = True

if not found_any:
    print("No string values found — checking for string 'nan' after string conversion...")
    for col in properties:
        mask = df[col].astype(str) == 'nan'
        bad = df[mask]
        if len(bad) > 0:
            print(f"Column '{col}' has {len(bad)} string-nan rows:")
            print(bad[['sequence', col]])
