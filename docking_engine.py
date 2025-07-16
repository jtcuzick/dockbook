import os
import json
import subprocess
import streamlit as st

def clean_ligand_file(path):
    with open(path, 'r') as file:
        lines = file.readlines()
    cleaned = [line for line in lines if not line.startswith("MODEL") and not line.startswith("ENDMDL")]
    with open(path, 'w') as file:
        file.writelines(cleaned)

def run_basic_docking(vina_path, receptor_path, ligand_folder, output_folder, log_folder, center, size):
    receptor_name = os.path.splitext(os.path.basename(receptor_path))[0]
    os.makedirs(output_folder, exist_ok=True)
    os.makedirs(log_folder, exist_ok=True)

    docking_scores = []
    docked_ligands = []
    
    # Reset scores list
    st.session_state["docking_scores"] = []

    for ligand_file in os.listdir(ligand_folder):
        if ligand_file.endswith(".pdbqt"):
            ligand_path = os.path.join(ligand_folder, ligand_file)
            clean_ligand_file(ligand_path)

            ligand_name = os.path.splitext(ligand_file)[0]
            receptor_name = os.path.splitext(os.path.basename(receptor_path))[0]
            out_file = os.path.join(output_folder, f"{ligand_name}_to_{receptor_name}.pdbqt")
            log_file = os.path.join(log_folder, f"{ligand_name}_to_{receptor_name}.txt")

            command = [
                vina_path,
                "--receptor", receptor_path,
                "--ligand", ligand_path,
                "--center_x", str(center[0]),
                "--center_y", str(center[1]),
                "--center_z", str(center[2]),
                "--size_x", str(size[0]),
                "--size_y", str(size[1]),
                "--size_z", str(size[2]),
                "--out", out_file
            ]

            result = subprocess.run(command, capture_output=True, text=True)

            # Save raw log
            with open(log_file, 'w') as log:
                log.write(result.stdout + result.stderr)

            # Extract score from stdout
            score = None
            for line in result.stdout.splitlines():
                if line.strip().startswith("1 "):  # mode 1: best pose
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            score = float(parts[1])
                        except ValueError:
                            score = None
                    break

            docking_scores.append((ligand_name, score))
            docked_ligands.append(out_file)

    # Save to session
    st.session_state["docking_scores"] = docking_scores

    # Save notebook-specific memory
    notebook_name = st.session_state.get("current_notebook", "Untitled")
    view_file = f"notebooks/{notebook_name}_view.json"

    memory_data = {
        "receptor": receptor_path,
        "ligands": docked_ligands,
        "scores": docking_scores
    }

    with open(view_file, 'w') as f:
        json.dump(memory_data, f)

