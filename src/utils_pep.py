import torch
import math
import torch.nn as nn

# Split peptide sequence into amino acid tokens
def split(seq):
    '''
    function: Split a peptide sequence into individual amino acid tokens.
    input: A peptide sequence string, e.g. "ACDEFG"
    output: A string with space between residues, e.g. "A C D E F G"
    '''
    seq = str(seq).strip().upper()
    return ' '.join(list(seq))

class GELU(nn.Module):
    def forward(self, x):
        return 0.5 * x * (1 + torch.tanh(math.sqrt(2 / math.pi) * (x + 0.044715 * torch.pow(x, 3))))

class PositionwiseFeedForward(nn.Module):
    def __init__(self, d_model, d_ff, dropout=0.1):
        super(PositionwiseFeedForward, self).__init__()
        self.w_1 = nn.Linear(d_model, d_ff)
        self.w_2 = nn.Linear(d_ff, d_model)
        self.dropout = nn.Dropout(dropout)
        self.activation = GELU()

    def forward(self, x):
        return self.w_2(self.dropout(self.activation(self.w_1(x))))

class LayerNorm(nn.Module):
    def __init__(self, features, eps=1e-6):
        super(LayerNorm, self).__init__()
        self.a_2 = nn.Parameter(torch.ones(features))
        self.b_2 = nn.Parameter(torch.zeros(features))
        self.eps = eps

    def forward(self, x):
        mean = x.mean(-1, keepdim=True)
        std = x.std(-1, keepdim=True)
        return self.a_2 * (x - mean) / (std + self.eps) + self.b_2

class SublayerConnection(nn.Module):
    def __init__(self, size, dropout):
        super(SublayerConnection, self).__init__()
        self.norm = LayerNorm(size)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x, sublayer):
        return x + self.dropout(sublayer(self.norm(x)))

# Sample peptide tokens from probabilistic distribution
def sample(msms):
    ret = []
    for msm in msms:
        ret.append(torch.multinomial(msm.exp(), 1).squeeze())
    return torch.stack(ret)

# Check validity of a peptide sequence (standard amino acids only)
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

def validity(sequences):
    loss = 0
    for seq in sequences:
        seq = str(seq).strip().upper()
        if len(seq) == 0 or not all(aa in AMINO_ACIDS for aa in seq):
            loss += 1
    return 1 - loss / len(sequences)
