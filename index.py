import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from rdkit import Chem
from rdkit.Chem import Draw
from utils.generate_peptides import generate_peptides
from maccs_keys.index import maccs_keys
from kcf.index import kcfs
from donor_acceptor.index import donor_accepter

def main(target, peptide_length, method):
    Draw.MolToImage(Chem.MolFromSmiles(target)).save("./output/target.png")

    print("Generating peptides...")
    peptides = generate_peptides(peptide_length)
    print(f"Generated {len(peptides)} peptides.")
    
    if method == 'maccs_keys':
        results = maccs_keys(target, peptides)
    elif method == 'kcfs':
        results = kcfs(target, peptides)
    elif method == 'donor_accepter':
        results = donor_accepter(target, peptides)
    else:
        raise ValueError("Invalid method. Choose either 'maccs_keys' or 'kcfs'.")
    
    # csvファイルに結果を書き込む
    with open("./output/results.csv", "w") as f:
        f.write("peptide,similarity\n")
        for peptide, similarity in results:
            f.write(f"{peptide},{similarity}\n")
    
    # 類似度の高いペプチドを描画
    high_score_peptides = [result[0] for result in results[:9]]
    rdkit_high_score_peptides = [Chem.MolFromFASTA(peptide) for peptide in high_score_peptides]
    legends = [f"{peptide} ({similarity:.3f})" for peptide, similarity in results[:9]]
    Draw.MolsToGridImage(rdkit_high_score_peptides, molsPerRow=3, legends=legends).save("./output/high_score_peptides.png")

    print("All processes are finished. Check the output folder for the results.")

if __name__ == "__main__":
    target = "CCN(C(=O)C1=C(C)ON=C1C1=C(Cl)C=CC=C1)C1=CC(Cl)=CC=C1"
    peptide_length = 3
    
    print("Choose the method:")
    print("1. maccs_keys")
    print("2. kcfs")
    print("3. donor_accepter")
    
    method_choice = input("Enter the number of the method: ")
    if method_choice == "1":
        method = "maccs_keys"
    elif method_choice == "2":
        method = "kcfs"
    elif method_choice == "3":
        method = "donor_accepter"
    else:
        print("Invalid choice. Defaulting to maccs_keys.")
        method = "maccs_keys"
    
    main(target, peptide_length, method)