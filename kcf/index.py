"""
MACCS Keysを使った類似度計算

@params
target: 比較対象の分子のSMILES表記
peptides: 比較するペプチドのリスト

@return
sorted_results: 類似度の高い順にソートされた結果。型はリスト。各要素はタプルで、(ペプチド, 類似度)の形式。
"""
import sys
sys.path.append("~/SAR")

from . import kcfconvoy as kcf
from rdkit import Chem
from .kcfconvoy import similarity
import time
import logging

logging.basicConfig(level=logging.INFO)

def kcfs(target, peptides):

    print('Converting target molecule to KCF vector...')

    kcf_target = kcf.KCFvec()
    kcf_target.input_smiles(target)
    kcf_target.convert_kcf_vec()

    print(f'Processing {len(peptides)} peptides...')
    
    results = []
    start_time = time.time()
    total_peptides = len(peptides)
    
    for i, peptide in enumerate(peptides):
        smiles_peptide = Chem.MolToSmiles(Chem.MolFromFASTA(peptide))
        kcf_peptide = kcf.KCFvec()
        kcf_peptide.input_smiles(smiles_peptide)
        kcf_peptide.convert_kcf_vec()

        # similarity_resultは3つの要素からなるタプルで、(重み付きタニモト係数, peptide中のtargetが持つ部分構造の比率、target中のpeptideが持つ部分構造の比率)の形式。
        similarity_result = similarity(kcf_target, kcf_peptide)
        results.append((peptide, similarity_result[0]))

        # ログ出力
        elapsed_time = time.time() - start_time
        progress = (i + 1) / total_peptides * 100
        logging.info(f'Processed {i + 1}/{total_peptides} peptides ({progress:.2f}%). Elapsed time: {elapsed_time:.2f} seconds.')

    sorted_results = sorted(results, key=lambda x: x[1], reverse=True)
    return sorted_results

if __name__ == "__main__":
    target = "CCC1=C(N(N=N1)C1=CC=CC=C1)C(=O)N(C)C1=CC=C(Cl)C=C1"
    peptides = ['HIP', 'IPH', 'HIF', 'IHF']
    print(kcfs(target, peptides))
    pass