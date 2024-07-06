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
    name = "Oximes"
    target = "CCN(C(=O)C1=C(C)ON=C1C1=C(Cl)C=CC=C1)C1=CC(Cl)=CC=C1"
    peptide_length = 3
    
    main(name, target, peptide_length)