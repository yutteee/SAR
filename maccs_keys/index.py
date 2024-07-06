"""
MACCS Keysを使った類似度計算

@params
target: 比較対象の分子のSMILES表記
peptides: 比較するペプチドのリスト

@return
sorted_results: 類似度の高い順にソートされた結果。型はリスト。各要素はタプルで、(ペプチド, 類似度)の形式。
"""

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

def maccs_keys(target, peptides):
    print('Converting target molecule to MACCS keys...')
    rdkit_target = Chem.MolFromSmiles(target)
    maccs_target = AllChem.GetMACCSKeysFingerprint(rdkit_target)
    
    print(f'Processing {len(peptides)} peptides...')
    maccs_peptides = []
    for peptide in peptides:
        rdkit_peptide = Chem.MolFromFASTA(peptide)
        maccs_peptide = AllChem.GetMACCSKeysFingerprint(rdkit_peptide)
        maccs_peptides.append(maccs_peptide)

    print('Calculating Tanimoto similarity...')
    maccs = DataStructs.BulkTanimotoSimilarity(maccs_target, maccs_peptides)
    results = dict(zip(peptides, maccs))
    sorted_results = sorted(results.items(), key=lambda x: x[1], reverse=True)
    return sorted_results

if __name__ == "__main__":
    target = "CCC1=C(N(N=N1)C1=CC=CC=C1)C(=O)N(C)C1=CC=C(Cl)C=C1"
    peptides = ['HIP', 'IPH', 'HIF', 'IHF']
    print(maccs_keys(target, peptides))
    pass