import pandas as pd
df = pd.read_csv('../data/peptide_conditional_norm.csv')
print(df.shape)
print(df[['sweet','bitter','umami','sour','salty','length','charge','mw','aromatic','gravy','hydrophobic_fraction']].dtypes)
for col in ['sweet','bitter','umami','sour','salty','length','charge','mw','aromatic','gravy','hydrophobic_fraction']:
    print(col, df[col].dtype)
    if df[col].dtype == object:
        print(df[col].unique()[:20])
        print(sorted(['sweet','bitter','umami','sour','salty','length','charge','mw','aromatic','gravy','hydrophobic_fraction']))
        print(list(df.columns))
