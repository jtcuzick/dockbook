from flask import Flask, request, jsonify
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
import tempfile
import os

app = Flask(__name__)

@app.route('/smiles-to-sdf', methods=['POST'])
def smiles_to_sdf():
    data = request.json
    smiles = data.get("smiles")
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return jsonify({"error": "Invalid SMILES"}), 400

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
    writer = Chem.SDWriter(tmp.name)
    writer.write(mol)
    writer.close()

    with open(tmp.name, "r") as f:
        sdf_data = f.read()
    os.unlink(tmp.name)

    return jsonify({"sdf": sdf_data})

@app.route('/sdf-to-pdbqt', methods=['POST'])
def sdf_to_pdbqt():
    data = request.json
    sdf_content = data.get("sdf")

    with tempfile.NamedTemporaryFile(delete=False, suffix=".sdf") as sdf_file:
        sdf_file.write(sdf_content.encode())
        sdf_path = sdf_file.name

    pdbqt_path = sdf_path.replace(".sdf", ".pdbqt")

    result = subprocess.run([
        "obabel", sdf_path, "-O", pdbqt_path,
        "--partialcharge", "eem", "--gen3d", "--addh"
    ], capture_output=True)

    if result.returncode != 0:
        return jsonify({"error": result.stderr.decode()}), 500

    with open(pdbqt_path, "r") as f:
        pdbqt_data = f.read()

    os.unlink(sdf_path)
    os.unlink(pdbqt_path)

    return jsonify({"pdbqt": pdbqt_data})

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
