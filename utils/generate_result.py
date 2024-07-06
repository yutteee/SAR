"""
結果をoutputディレクトリに保存する関数
"""
import os
from rdkit import Chem
from rdkit.Chem import Draw

def generate_result(target, results, output_dir="./output"):
    os.makedirs(output_dir, exist_ok=True)

    # targetの構造を描画
    Draw.MolToImage(Chem.MolFromSmiles(target)).save(f'{output_dir}/target.png')

    # csvファイルに結果を書き込む
    with open(f'{output_dir}/result.csv', "w") as f:
        f.write("peptide,similarity\n")
        for peptide, similarity in results:
            f.write(f"{peptide},{similarity}\n")
    
    # 類似度の高いペプチドを描画
    high_score_peptides = [result[0] for result in results[:9]]
    rdkit_high_score_peptides = [Chem.MolFromFASTA(peptide) for peptide in high_score_peptides]
    legends = [f"{peptide} ({similarity:.3f})" for peptide, similarity in results[:9]]
    Draw.MolsToGridImage(rdkit_high_score_peptides, molsPerRow=3, legends=legends).save(f'{output_dir}/high_score_peptides.png')

if __name__ == "__main__":
    target = "CCN(C(=O)C1=C(C)ON=C1C1=C(Cl)C=CC=C1)C1=CC(Cl)=CC=C1"
    results = [('HIP', 0.8), ('IPH', 0.7), ('HIF', 0.6), ('IHF', 0.5)]
    generate_result(target, results, "../output")
    pass