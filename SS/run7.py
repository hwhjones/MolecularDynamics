import streamlit as st
import plotly.graph_objects as go
from Bio import PDB
import numpy as np

st.title("PDB 3D Viewer")

# Path to your local PDB file
pdb_file_path = "converted.pdb"

# Parse the PDB file
parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure("protein", pdb_file_path)

# Get the number of models
num_models = len(structure)
st.write(f"Number of models: {num_models}")

# Function to extract coordinates from a model
def get_model_coordinates(model):
    coords = []
    for chain in model:
        for residue in chain:
            for atom in residue:
                coords.append(atom.coord)
    return np.array(coords)

# Frame selection slider
frame = st.slider("Select frame", 0, num_models-1, 0)

# Get coordinates for the selected model
coords = get_model_coordinates(structure[frame])

# Create 3D scatter plot
fig = go.Figure(data=[go.Scatter3d(
    x=coords[:, 0],
    y=coords[:, 1],
    z=coords[:, 2],
    mode='markers',
    marker=dict(
        size=2,
        color=coords[:, 2],  # Color by z-coordinate
        colorscale='Viridis',
        opacity=0.8
    )
)])

# Update layout for better visibility
fig.update_layout(
    scene=dict(
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z',
        aspectmode='data'
    ),
    width=700,
    height=700,
    title=f"Model {frame + 1}"
)

# Display the plot
st.plotly_chart(fig)

# Button to cycle through frames
if st.button("Next Frame"):
    frame = (frame + 1) % num_models
    st.session_state.frame = frame
    st.experimental_rerun()

# Initialize session state for frame if not exists
if 'frame' not in st.session_state:
    st.session_state.frame = 0

# Auto-play option
auto_play = st.checkbox("Auto-play")

if auto_play:
    import time
    frame = st.session_state.frame
    frame = (frame + 1) % num_models
    st.session_state.frame = frame
    time.sleep(0.5)  # Adjust this value to control animation speed
    st.experimental_rerun()