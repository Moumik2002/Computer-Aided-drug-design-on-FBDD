
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

# Example molecule: Aspirin
smiles = "CC(=O)Oc1ccccc1C(=O)O"

mol = Chem.MolFromSmiles(smiles)

if mol is None:
    print("Invalid SMILES")
else:
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)

    violations = sum([
        mw > 500,
        logp > 5,
        hbd > 5,
        hba > 10
    ])

    print("Molecular weight:", round(mw, 2))
    print("LogP:", round(logp, 2))
    print("H-bond donors:", hbd)
    print("H-bond acceptors:", hba)
    print("Lipinski violations:", violations)

    if violations <= 1:
        print("Passes Lipinski Rule of Five")
    else:
        print("Fails Lipinski Rule of Five")

