import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from utils.generate_peptides import generate_peptides
from utils.generate_result import generate_result
from maccs_keys.index import maccs_keys
from kcf.index import kcfs
from donor_acceptor.index import donor_accepter

def main(name, target, peptide_length):
    peptides = generate_peptides(peptide_length)

    # MACCS Keys
    results = maccs_keys(target, peptides)
    generate_result(target, results, f'./output/{name}/maccs_keys')

    # KCFS
    results = kcfs(target, peptides)
    generate_result(target, results, f'./output/{name}/kcf')

    # Donor-Acceptor
    results = donor_accepter(target, peptides)
    generate_result(target, results, f'./output/{name}/donor_acceptor')

    print("All processes are finished. Check the output folder for the results.")

if __name__ == "__main__":
    peptide_length = 3

    target_pair = [
        {
            "name": "2-Aryl-3-aminomethylquinolines",
            "target": "COC1=CC2=C(C=C1)C=C(CNCCC1=CC=C(Br)C=C1)C(=N2)C1=CSC=C1",
        },
        {
            "name": "Sulfonamide",
            "target": "COC1=CC=C(C=C1OC)S(=O)(=O)N1CCCN(CCC1)S(=O)(=O)C1=CC(OC)=C(OC)C=C1"
        },
        {
            "name": "Azepine",
            "target": "FC1=CC=CC(F)=C1C1=C(CN2[C@H](ON=C2C2=CC=CO2)C2=CC(=CC(=C2Cl)C(F)(F)F)C(F)(F)F)C=NC=C1"
        },
        {
            "name": "4-Phenylpyridine",
            "target": "ClC1=CC(OC2=C(C=CC=C2)[N]2=CN=CC2C2=C(Cl)C=CC=C2)=C(Cl)C=C1"
        },
    ]
    
    for target in target_pair:
        main(target["name"], target["target"], peptide_length)
        pass