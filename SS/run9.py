import streamlit as st
import py3Dmol
from Bio import PDB
import io

# Create a Streamlit app
st.title("PDB Structure Viewer")

# Specify the path to your local PDB file
pdb_file_path = 'converted.pdb'  # Change this to your local PDB file path

# Read the PDB file
with open(pdb_file_path, 'r') as pdb_file:
    pdb_string = pdb_file.read()

# Parse the PDB file
parser = PDB.PDBParser()
structure = parser.get_structure("protein", io.StringIO(pdb_string))

# Create a 3Dmol viewer
viewer = py3Dmol.view(width=600, height=400)
viewer.addModel(pdb_string, 'pdb')
viewer.setStyle({'cartoon': {'color': 'spectrum'}})
viewer.zoomTo()
viewer.render()

# Count the number of models (frames)
num_models = len(list(structure.get_models()))

if num_models > 1:
    # Add a frame slider if there are multiple models
    frame = st.slider('Frame', 0, num_models - 1, 0)
    viewer.setFrame(frame)
    st.write(f"Current frame: {frame + 1} of {num_models}")

# Get the viewer HTML
viewer_html = viewer.render()

# Display the viewer using a custom component
st.components.v1.html(f"""
    <div style="width:620px;height:420px;">
        {viewer_html}
    </div>
""", height=420, width=620)

# Display additional information
st.write(f"Number of atoms: {len(list(structure.get_atoms()))}")
st.write(f"Number of residues: {len(list(structure.get_residues()))}")
st.write(f"Number of chains: {len(list(structure.get_chains()))}")