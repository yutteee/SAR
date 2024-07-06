"""
ドナーアクセプターペアフィンガープリントを使った類似度計算

@params
target: 比較対象の分子のSMILES表記
peptides: 比較するペプチドのリスト

@return
sorted_results: 類似度の高い順にソートされた結果。型はリスト。各要素はタプルで、(ペプチド, 類似度)の形式。
"""

from rdkit import Chem, DataStructs
from rdkit.Chem.AtomPairs import Sheridan

def donor_accepter(target, peptides):
    print('Converting target molecule to donor_acceptor_pair...')
    rdkit_target = Chem.MolFromSmiles(target)
    da_target = Sheridan.GetBPFingerprint(rdkit_target)
    
    print(f'Processing {len(peptides)} peptides...')
    da_peptides = []
    for peptide in peptides:
        rdkit_peptide = Chem.MolFromFASTA(peptide)
        da_peptide = Sheridan.GetBPFingerprint(rdkit_peptide)
        da_peptides.append(da_peptide)

    print('Calculating Tanimoto similarity...')
    da = DataStructs.BulkTanimotoSimilarity(da_target, da_peptides)
    results = dict(zip(peptides, da))
    sorted_results = sorted(results.items(), key=lambda x: x[1], reverse=True)
    return sorted_results

if __name__ == "__main__":
    target = "CCC1=C(N(N=N1)C1=CC=CC=C1)C(=O)N(C)C1=CC=C(Cl)C=C1"
    peptides = ['HIP', 'IPH', 'HIF', 'IHF']
    print(donor_accepter(target, peptides))
    pass