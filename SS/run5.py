import streamlit as st
import py3Dmol
import streamlit.components.v1 as components
from stmol import showmol

st.title("PDB Viewer")

# Path to your local PDB file
pdb_file_path = "converted.pdb"

# Read the PDB file
with open(pdb_file_path, 'r') as f:
    pdb_content = f.read()

# Create a py3Dmol view
view = py3Dmol.view(width=700, height=500)
view.addModel(pdb_content, "pdb")
view.setStyle({'cartoon': {'color': 'spectrum'}})
view.zoomTo()
showmol(view, height=400, width=800)