# DockBook

DockBook is a Streamlit-based web application for structure-based virtual screening. It integrates RDKit for cheminformatics, Open Babel for format conversion and 3D preparation, and py3Dmol for interactive molecular visualization and docking.

## Features

 Run Docking Workflows: Prepare ligands and receptors for docking simulations.

 Format Conversion: Seamless ligand file conversion (SMILES, SDF, PDBQT, etc.) via Open Babel.

 3D Visualization: Interactive molecule rendering powered by py3Dmol.

 Batch Screening Support: Process multiple ligands in parallel and export results for downstream analysis in DataWarrior.

## How It Works

- Step 1: Upload receptor and ligand files.  
- Step 2: Ligands are converted into 3D with Open Babel.  
- Step 3: Docking workflow is run (AutoDock Vina or future engines).  
- Step 4: Results are visualized interactively with py3Dmol.
- Step 5: Export results in .txt file for analysis in DataWarrior.

## Access

You can try DockBook directly in your browser here:  

 [Launch DockBook on Streamlit Cloud](https://dockbook-v1.streamlit.app/)

## Installation (Developers)

Clone the repository and install dependencies:

```bash
git clone https://github.com/<your-username>/dockbook.git
cd dockbook
pip install -r requirements.txt



