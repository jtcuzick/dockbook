import streamlit as st
import os
import json
from docking_engine import run_basic_docking
from viewer import render_docked_molecules
import streamlit.components.v1 as components
import tempfile
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolfiles
import subprocess
import threading

# Ensure a notebook is selected
if "current_notebook" not in st.session_state or not st.session_state["current_notebook"]:
    st.error("❌ No notebook selected. Please go back to the home screen and open a notebook.")
    st.stop()

# Prepare notebook name and paths
notebook_name = st.session_state["current_notebook"]
base_path = os.path.join("notebooks", notebook_name)
notebook_file = os.path.join("notebooks", f"{notebook_name}.json")
view_file = os.path.join(base_path, "view.json")

# Create necessary directories for this notebook
receptor_path = os.path.join(base_path, "receptor.pdbqt")
ligand_folder = os.path.join(base_path, "ligands")
output_folder = os.path.join(base_path, "docked")
log_folder = os.path.join(base_path, "logs")

for folder in [base_path, ligand_folder, output_folder, log_folder]:
    os.makedirs(folder, exist_ok=True)

# Load saved notebook config (UI fields)
if os.path.exists(notebook_file):
    with open(notebook_file, 'r') as f:
        saved_data = json.load(f)
        for key, value in saved_data.items():
            st.session_state[key] = value

col1, col2, col3 = st.columns([1, 6, 1])
with col1:
    if st.button("← Back", key="back_button"):
        st.session_state.current_notebook = None
        st.switch_page("app.py")

with col3:
    if os.path.exists(view_file):
        if st.button("Next →", key="forward_button"):
            st.switch_page("pages/ViewDocking.py")


# User Inputs

st.title(f" NOTEBOOK: {notebook_name}")
st.markdown("Set up your docking run by completing the following setup mechanisms.")

st.markdown("""
<div style='color: #4a90e2; font-size: 1.75em; font-weight: 700;'>
    1) Ligand Builder
</div>
""", unsafe_allow_html=True)

st.markdown("""
Draw your desired ligand **below**, copy the respective SMILES string, and a 3D structure will be generated and converted to **PDBQT** format for docking.
""")


# Embedding Kekule.js editor with SMILES export logic
components.html("""
<!DOCTYPE html>
<html>
<head>
  <script src="https://unpkg.com/kekule@1.0.2/dist/kekule.min.js"></script>
  <link rel="stylesheet" href="https://unpkg.com/kekule@1.0.2/dist/themes/default/kekule.css">
  <style>
    #chemEditor {
      width: 100%;
      height: 500px;
      border: 1px solid #ccc;
      margin-bottom: 1em;
    }
    #output {
      font-family: monospace;
      padding: 0.5em;
      background-color: #f0f0f0;
    }
  </style>
</head>
<body>
  <div id="chemEditor"></div>
  <button onclick="getSmiles()">Get SMILES</button>
  <p id="output">SMILES will appear here...</p>
  <script>
    var composer = new Kekule.Editor.Composer(document.getElementById('chemEditor'));
    composer.setPredefinedSetting('fullFunc');

    function getSmiles() {
      const mol = composer.getChemObj();
      if (mol) {
        const smiles = Kekule.IO.saveFormatData(mol, 'smi');
        document.getElementById("output").innerText = smiles;
      } else {
        document.getElementById("output").innerText = "No molecule drawn.";
      }
    }
  </script>
</body>
</html>
""", height=600)


# Input SMILES
smiles = st.text_input("Paste SMILES string here")

# Function to generate RDKit 3D molecule
def smiles_to_3d_mol(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    return mol

# Function to save RDKit molecule as PDB
def save_mol_as_pdb(mol, output_path: str):
    writer = rdmolfiles.PDBWriter(output_path)
    writer.write(mol)
    writer.close()

# Convert PDB to PDBQT using Open Babel
def convert_pdb_to_pdbqt(pdb_path: str, pdbqt_path: str):
    command = ["obabel", pdb_path, "-O", pdbqt_path, "--gen3d"]
    result = subprocess.run(command, capture_output=True, text=True)
    return result.returncode == 0, result.stdout + result.stderr

if st.button("Generate PDBQT"):
    if not smiles:
        st.error("❌ Please provide a valid SMILES string.")
    else:
        mol = smiles_to_3d_mol(smiles)
        if mol is None:
            st.error("❌ Invalid SMILES string.")
        else:
            with tempfile.TemporaryDirectory() as tmpdir:
                pdb_path = os.path.join(tmpdir, "ligand.pdb")
                pdbqt_path = os.path.join(tmpdir, "ligand.pdbqt")

                save_mol_as_pdb(mol, pdb_path)
                success, output = convert_pdb_to_pdbqt(pdb_path, pdbqt_path)

                if success and os.path.exists(pdbqt_path):
                    with open(pdbqt_path, "r") as f:
                        pdbqt_contents = f.read()
                    st.success("✅ PDBQT generated!")
                    st.download_button("Download PDBQT", pdbqt_contents, file_name="ligand.pdbqt")
                else:
                    st.error("❌ Conversion to PDBQT failed.")
                    st.text(output)
                    
st.caption("⚠️ IMPORTANT: Due to browser restrictions, the downloaded file will be saved to your system’s default download folder. Please move it manually to your desired directory for later use in the docking workflow.")

st.markdown("""
<div style='color: #4a90e2; font-size: 1.75em; font-weight: 700;'>
    2) Establish File Pathways
</div>
""", unsafe_allow_html=True)                    


st.markdown("""
Enter your file pathways for the **AutoDock Vina Executable** (see licensing details below), **receptor file** (usually retrieved from Protein Data Base), and **ligand folder** (gather all saved PDBQT files from the previous step in one file). 
""")

vina_path = st.text_input(" Path to AutoDock Vina Executable", value=st.session_state.get("vina_path", ""))

st.caption("""
Working in accordance with the Apache License. For more information please visit: http://www.apache.org/licenses/
""")

receptor_path = st.text_input(" Path to Receptor (.pdbqt)", value=st.session_state.get("receptor_path", ""))
ligand_folder = st.text_input(" Path to Ligand Folder", value=st.session_state.get("ligand_folder", ""))
st.session_state["ligand_input_folder"] = ligand_folder

st.markdown("""
<div style='color: #4a90e2; font-size: 1.75em; font-weight: 700;'>
    3) Center Coordinates and Docking Box
</div>
""", unsafe_allow_html=True)

st.markdown("""
Enter these measurements to the nearest hundredths place.
""")

center_x = st.number_input("Center X (Å)", value=float(st.session_state.get("center_x", 0.0)), step=0.1)
center_y = st.number_input("Center Y (Å)", value=float(st.session_state.get("center_y", 0.0)), step=0.1)
center_z = st.number_input("Center Z (Å)", value=float(st.session_state.get("center_z", 0.0)), step=0.1)

size_x = st.number_input("Docking Box - Size X (Å)", value=float(st.session_state.get("size_x", 0.0)), step=0.1)
size_y = st.number_input("Docking Box - Size Y (Å)", value=float(st.session_state.get("size_y", 0.0)), step=0.1)
size_z = st.number_input("Docking Box - Size Z (Å)", value=float(st.session_state.get("size_z", 0.0)), step=0.1)


# Save Notebook State
to_save = {
    "vina_path": vina_path,
    "receptor_path": receptor_path,
    "ligand_folder": ligand_folder,
    "center_x": center_x,
    "center_y": center_y,
    "center_z": center_z,
    "size_x": size_x,
    "size_y": size_y,
    "size_z": size_z
}

for key, value in to_save.items():
    st.session_state[key] = value

with open(notebook_file, 'w') as f:
    json.dump(to_save, f)

# Docking Run
output_folder = os.path.join("results/docked_structures", notebook_name)
log_folder = os.path.join("results/docking_logs", notebook_name)
os.makedirs(output_folder, exist_ok=True)
os.makedirs(log_folder, exist_ok=True)

def get_docking_score(ligand_name):
    # Extract original ligand name
    basename = os.path.basename(ligand_name)
    ligand_prefix = basename.split("_to_")[0]  # Get 'Ligand_1' from 'Ligand_1_to_Receptor.pdbqt'
    
    # Use notebook-specific log folder
    notebook_name = st.session_state.get("current_notebook", "Untitled")
    log_file = os.path.join("results/docking_logs", notebook_name, f"{ligand_prefix}.txt")

    if not os.path.exists(log_file):
        return "N/A"

    with open(log_file, 'r') as f:
        for line in f:
            if "REMARK VINA RESULT:" in line:
                try:
                    return float(line.split()[3])
                except:
                    return "N/A"
    return "N/A"


if st.button("🚀 RUN"):
    if not os.path.isfile(vina_path):
        st.error("❌ Invalid Vina executable path.")
    elif not os.path.isfile(receptor_path):
        st.error("❌ Invalid receptor file path.")
    elif not os.path.isdir(ligand_folder):
        st.error("❌ Invalid ligand folder path.")
    else:
        
        # Save empty state and move to visualization
        st.session_state["viewer_receptor"] = receptor_path
        st.session_state["viewer_ligands"] = []
        st.session_state["docking_scores"] = []

        view_file = os.path.join(base_path, "view.json")
        with open(view_file, 'w') as f:
            json.dump({
                "receptor": receptor_path,
                "ligands": [],
                "scores": []
            }, f)

        # Start docking in background
        threading.Thread(
            target=run_basic_docking,
            args=(
                vina_path,
                receptor_path,
                ligand_folder,
                output_folder,
                log_folder,
                (center_x, center_y, center_z),
                (size_x, size_y, size_z)
            ),
            daemon=True
        ).start()

        # Navigate to visualization
        st.switch_page("pages/ViewDocking.py")


