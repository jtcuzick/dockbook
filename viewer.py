import py3Dmol
import streamlit.components.v1 as components

def render_docked_molecules(receptor_file, ligand_files):
    width, height = 800, 500
    view = py3Dmol.view(width=width, height=height)

    # Load and style receptor
    with open(receptor_file, 'r') as f:
        receptor_data = f.read()
        view.addModel(receptor_data, "pdbqt")
        view.setStyle({'model': 0}, {'stick': {}})

    # Load and style ligands
    for i, ligand_file in enumerate(ligand_files):
        with open(ligand_file, 'r') as f:
            ligand_data = f.read()
            view.addModel(ligand_data, "pdbqt")
            view.setStyle({'model': i + 1}, {'stick': {'color': 'green'}})

    view.zoomTo()

    # Render the full HTML view
    full_html = view._make_html()
    
    # Just wrap the full HTML in a centered box
    components.html(
        f"""
        <div style="display: flex; justify-content: center;">
            <div style="border: 2px solid black; border-radius: 10px;">
                {full_html}
            </div>
        </div>
        """,
        height=height + 60,
        scrolling=False
    )




