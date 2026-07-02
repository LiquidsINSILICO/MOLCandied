import torch
import pandas as pd
from build_vocab import WordVocab
from utils_pep import split

SMI_MAX_SIZE = 200

with open("../data/peptides_corpus.txt", "r") as f:
    pep_vocab = WordVocab(f, max_size=SMI_MAX_SIZE, min_freq=1)
print("Vocab size:", len(pep_vocab))

from torch.utils.data import Dataset

class CustomTargetDataset(Dataset):
    def __init__(self, df, vocab, properties_list):
        self.df = df.reset_index(drop=True)
        self.vocab = vocab
        self.props = properties_list

    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        row = self.df.iloc[idx]
        pep = row['sequence']
        seq = self.vocab.to_seq(split(pep), seq_len=SMI_MAX_SIZE, with_eos=True, with_sos=True)
        prop_vals = [row[p] for p in self.props]
        return {
            "peptide_rep": torch.tensor(seq, dtype=torch.long),
            "properties":  torch.tensor(prop_vals, dtype=torch.float32),
            "sequence":    pep
        }

df = pd.read_csv('../data/peptide_conditional_norm.csv')
train_df = df[df['source'] == 'train'].reset_index(drop=True)

properties = ['sweet','bitter','umami','sour','salty','length','charge',
              'mw','aromatic','gravy','hydrophobic_fraction']

dataset = CustomTargetDataset(train_df, pep_vocab, properties_list=properties)

print(f"Testing {len(dataset)} rows...")
for i in range(len(dataset)):
    try:
        item = dataset[i]
        if item['properties'].dtype != torch.float32:
            print(f"Row {i} has wrong properties dtype: {item['properties'].dtype}")
        if item['peptide_rep'].dtype != torch.long:
            print(f"Row {i} has wrong peptide_rep dtype: {item['peptide_rep'].dtype}")
    except Exception as e:
        print(f"Row {i} FAILED: {e}")
        print(train_df.iloc[i])
        break
else:
    print("All rows processed successfully — no errors found")
