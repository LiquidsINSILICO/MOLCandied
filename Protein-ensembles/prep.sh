#!/bin/bash
for i in {0..98}; do
        input="apo_rep100_c${i}.pdb"
        output="apo_rep100_c${i}_clean.pdb"
#================================================================================
# check if input file exist
#================================================================================
        if [ -f "$input" ]; then
                grep -v POPC "$input" > "$output"
                echo "processed: $input --> $output"
        else
                echo "Warning: $input NOT FOUND"
        fi
done
#
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
#================================================================================
#""Map the binding pocket residues from 9UTB pdb to simulated structure pdb
#   Use the python code created by Anupama""
#================================================================================
python map_pocket-fixed.py

# Define active site residue
echo "" > protein.act-pass
echo "42 46 78 117 118 140 252 253 277 359" >> protein.act-pass
#================================================================================
# ""You might come accross error for the residue name and atom name""
# ""Few of them for which haddock3-restraints threw error were modified to standard""
#================================================================================
sed -e 's/HID/HIS/g' -e 's/HSD/HIS/g' -e 's/CD  ILE /CD1 ILE /g' -e 's/OT1/OXT/g' -e 's/OT2/O  /g' Model1.pdb > Model1_fixed.pdb
#================================================================================
# ""Protein passive residue list were defined around active residues that are within radii of 7.5""
#================================================================================
haddock3-restraints passive_from_active Model1_fixed.pdb 42,46,78,117,118,140,252,253,277,359 -r 7.5 > protein_passive.list
#================================================================================
# ""Follow the similar rule for the protein ensemble structure as well""
#================================================================================
sed -e 's/HES/HIS/g' -e 's/HID/HIS/g' -e 's/HSD/HIS/g' -e 's/CD  ILE /CD1 ILE /g' -e 's/OT1/OXT/g' -e 's/OT2/O  /g' receptor_ensemble-V2.pdb > receptor_ensemble-V2_fixed.pdb
#================================================================================
#
#================================================================================
