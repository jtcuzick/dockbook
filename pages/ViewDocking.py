import streamlit as st
from viewer import render_docked_molecules
import os
import json
import time

# Set up
st.set_page_config(page_title="Docking Visualization", layout="wide")

# Load current notebook name
notebook_name = st.session_state.get("current_notebook", "Untitled")
output_dir = f"results/docked_structures/{notebook_name}"  # notebook-specific docking folder
view_file = f"notebooks/{notebook_name}_view.json"

# Load receptor and docking scores from memory (if available)
if os.path.exists(view_file):
    with open(view_file, 'r') as f:
        memory = json.load(f)
        st.session_state["viewer_receptor"] = memory.get("receptor", "")
        st.session_state["docking_scores"] = memory.get("scores", [])

# Ligand Polling

ligand_folder = st.session_state.get("ligand_input_folder")
if ligand_folder and os.path.exists(ligand_folder):
    expected_ligands = len([f for f in os.listdir(ligand_folder) if f.endswith(".pdbqt")])
    st.session_state["expected_ligand_count"] = expected_ligands

previous_ligand_count = len(st.session_state.get("viewer_ligands", []))
current_ligands = []
 
if os.path.exists(output_dir):
    current_ligands = sorted([
        os.path.join(output_dir, f)
        for f in os.listdir(output_dir)
        if f.endswith(".pdbqt")
    ])

# Update state if changed
if current_ligands != st.session_state.get("viewer_ligands", []):
    st.session_state["viewer_ligands"] = current_ligands
    st.rerun()

ligands = st.session_state["viewer_ligands"]

# Save current ligands into state
st.session_state["viewer_ligands"] = ligands

# Navigation Bar
col4, col5, col6 = st.columns([1, 6, 1])
with col4:
    if st.button("← Back", key="back_button"):
        st.switch_page("pages/DockingNotebook.py")
with col6:
    if st.button("Next →", key="forward_button"):
        st.switch_page("pages/Analysis.py")

# Visualization Section
st.markdown("""
    <div style='display: flex; flex-direction: column; align-items: center; justify-content: center; text-align: center;'>
        <h1 style='margin-bottom: 1rem;'>Docking Visualization</h1>
    </div>
""", unsafe_allow_html=True)

receptor = st.session_state.get("viewer_receptor")

if receptor:
    expected_count = st.session_state.get("expected_ligand_count")

    if ligands:
        if expected_count and len(ligands) >= expected_count:
            st.success(f"✅ Docking complete! {len(ligands)} ligands visualized.")
            
            ligand_names = [os.path.basename(l) for l in ligands]
            display_options = ["All Ligands"] + ligand_names

            selected = st.selectbox("Select ligand(s) to visualize:", display_options)

            if selected == "All Ligands":
                render_docked_molecules(receptor, ligands)
            else:
                selected_index = ligand_names.index(selected)
                render_docked_molecules(receptor, [ligands[selected_index]])
                
        else:
            st.info(f"⏳ Docking in progress... {len(ligands)} ligand(s) docked so far.")
            render_docked_molecules(receptor, ligands)  # show interim results
            time.sleep(5)
            st.rerun()
            
    else:
        st.info("⏳ Waiting for first ligand to dock...")
        time.sleep(5)
        st.rerun()
    
else:
    st.warning("🛑 No receptor loaded. Please start a docking run first.")







