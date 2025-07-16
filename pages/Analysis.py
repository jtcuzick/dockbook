import streamlit as st
import os
import subprocess

# Constants
notebook = st.session_state.get("current_notebook", "Untitled")
docking_dir = f"results/docked_structures/{notebook}"
results_file = f"results/{notebook}_vs_results.txt"

# Functions
def extract_smiles(input_file, output_file):
    command = f'/Users/jakecuzick/miniconda3/bin/obabel -ipdbqt {input_file} -osmi -O {output_file}'
    subprocess.run(command, shell=True, check=True)

def extract_docking_score(pdbqt_file):
    try:
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.strip().startswith("REMARK VINA RESULT:"):
                    parts = line.split()
                    return parts[3]
    except:
        return "N/A"

def generate_vs_results(directory, output_path):
    with open(output_path, 'w') as f:
        f.write("Title\tDocking Score\tSMILES\n")
        for file in sorted(os.listdir(directory)):
            if file.endswith(".pdbqt"):
                pdbqt_path = os.path.join(directory, file)
                smiles_path = os.path.splitext(pdbqt_path)[0] + ".smi"

                try:
                    extract_smiles(pdbqt_path, smiles_path)
                    with open(smiles_path, 'r') as smi:
                        smiles = smi.readline().strip().split()[0]
                except:
                    smiles = "N/A"

                score = extract_docking_score(pdbqt_path)
                title = os.path.splitext(file)[0]
                f.write(f"{title}\t{score}\t{smiles}\n")

# UI
st.title(" Docking Analysis")
st.markdown("Generate and export a `.txt` summary with docking scores and SMILES for use in **DataWarrior**.")

if st.button(" Generate Summary File"):
    try:
        generate_vs_results(docking_dir, results_file)
        with open(results_file, "r") as f:
            st.download_button("Download Results File", f.read(), file_name=os.path.basename(results_file))
    except Exception as e:
        st.error(f"❌ Failed to generate results file: {e}")


col1, col2 = st.columns([1, 9])
with col1:
    st.markdown("<div style='margin-top: 1em;'></div>", unsafe_allow_html=True)
    if st.button("← Back"):
        st.switch_page("pages/ViewDocking.py")

