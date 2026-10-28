from Bio.PDB import PDBParser
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import seq1

# 1. Load Structures
parser = PDBParser(QUIET=True)
ref_struct = parser.get_structure("9UTB", "9UTB.pdb")
model_struct = parser.get_structure("Model1", "Model1.pdb")

# Select Chain A (Change 'A' if your binding pocket is in a different chain)
ref_chain = ref_struct[0]['A']
model_chain = model_struct[0]['A']

# 2. Extract Sequences and Residue Lists
def get_seq_and_res_list(chain):
    seq = []
    res_list = []
    for res in chain:
        if res.has_id('CA'): # Only standard amino acids
            try:
                seq.append(seq1(res.get_resname()))
                res_list.append(res)
            except KeyError:
                continue # Skip non-standard residues
    return "".join(seq), res_list

ref_seq, ref_res_list = get_seq_and_res_list(ref_chain)
model_seq, model_res_list = get_seq_and_res_list(model_chain)

# 3. Perform Sequence Alignment using the NEW PairwiseAligner
aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -2
aligner.extend_gap_score = -0.5

alignments = aligner.align(ref_seq, model_seq)
best_alignment = alignments[0] # Take the best alignment

# 4. Build the Mapping Dictionary
# The alignment object maps indices between the two sequences
mapping = {}
ref_idx = 0
model_idx = 0

# Iterate through the aligned sequences
for i, (ref_char, model_char) in enumerate(zip(best_alignment[0], best_alignment[1])):
    if ref_char != '-':
        if model_char != '-':
            # Match: Map current ref residue to current model residue
            if ref_idx < len(ref_res_list) and model_idx < len(model_res_list):
                ref_res = ref_res_list[ref_idx]
                model_res = model_res_list[model_idx]
                # Store mapping: (Chain, ResNum) -> (Chain, ResNum)
                mapping[(ref_chain.id, ref_res.get_id()[1])] = (model_chain.id, model_res.get_id()[1])
            model_idx += 1
        # else: Gap in model (residue deleted in simulation)
        ref_idx += 1
    elif model_char != '-':
        # Gap in reference (insertion in simulation)
        model_idx += 1

# 5. Map Your Binding Pocket
binding_pocket_9UTB = [67, 71, 103, 142, 143, 165, 277, 278, 302, 384]   # change as you require

print(f"{'9UTB Residue':<15} | {'Model1 Residue':<15} | {'Status'}")
print("-" * 45)

for res_id in binding_pocket_9UTB:
    key = ('A', res_id) # Assuming Chain A
    if key in mapping:
        mapped_chain, mapped_resi = mapping[key]
        print(f"{res_id:<15} | {mapped_resi:<15} | Mapped (Chain {mapped_chain})")
    else:
        print(f"{res_id:<15} | {'Not Found':<15} | Gap or Missing Residue")   
