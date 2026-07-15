#!/bin/bash
for f in apo_rep100_c*_clean.pdb; do
    # 1. Extract 1-470 -> Chain A
    partA=$(pdb_selres -1:808 "$f" | pdb_chain -A)
    
    # 2. Extract 809-1291 -> Chain B
    partB=$(pdb_selres -809:1611 "$f" | pdb_chain -B)
    
    # 3. Merge and Tidy (pdb_tidy adds TER between chains automatically)
    echo "$partA" | cat - <(echo "$partB") | pdb_tidy -strict > "fixed_$f"
done

# Create ensemble
pdb_mkensemble fixed_apo_rep100_c*_clean.pdb | pdb_tidy -strict > receptor_ensemble-V2.pdb
rm fixed_apo_rep100_c*_clean.pdb   
