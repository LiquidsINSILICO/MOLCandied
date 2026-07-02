import torch
import pandas as pd
from build_vocab import WordVocab
from utils_pep import split
from torch.utils.data import Dataset

SMI_MAX_SIZE = 200

with open("../data/peptides_corpus.txt", "r") as f:
    pep_vocab = WordVocab(f, max_size=SMI_MAX_SIZE, min_freq=1)
print("Vocab size:", len(pep_vocab))

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
test_df  = df[df['source'] == 'test'].reset_index(drop=True)

properties = ['sweet','bitter','umami','sour','salty','length','charge',
              'mw','aromatic','gravy','hydrophobic_fraction']

for name, split_df in [('train', train_df), ('test', test_df)]:
    dataset = CustomTargetDataset(split_df, pep_vocab, properties_list=properties)
    print(f"\nTesting {name}: {len(dataset)} rows...")
    for i in range(len(dataset)):
        try:
            item = dataset[i]
        except Exception as e:
            print(f"{name} row {i} FAILED: {e}")
            print(split_df.iloc[i])
            break
    else:
        print(f"All {name} rows processed successfully")

from torch.utils.data import DataLoader

for name, split_df in [('train', train_df), ('test', test_df)]:
    dataset = CustomTargetDataset(split_df, pep_vocab, properties_list=properties)
    loader = DataLoader(dataset, batch_size=128, shuffle=True)
    print(f"\nTesting {name} DataLoader batches...")
    for i, batch in enumerate(loader):
        pass  # just iterate, let it crash if it's going to
    print(f"{name}: all {i+1} batches collated successfully")
