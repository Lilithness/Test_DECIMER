import sys, pathlib
from DECIMER.decimer import predict_SMILES       # <— NEW
from rdkit import Chem

img = pathlib.Path(sys.argv[1])
smiles = predict_SMILES(str(img))                # <— NEW
mol = Chem.MolFromSmiles(smiles)

print("SMILES:", smiles)
print("Valid:", mol is not None)
