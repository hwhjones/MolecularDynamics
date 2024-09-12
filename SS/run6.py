import streamlit as st
import nglview as nv
import MDAnalysis as mda
import tempfile

st.title("PDB Frame Viewer")

# Path to your local PDB file
pdb_file_path = "converted.pdb"

# Load the PDB file
u = mda.Universe(pdb_file_path)

# Get the number of frames
num_frames = len(u.trajectory)

st.write(f"Number of frames: {num_frames}")

# Create an nglview widget
view = nv.show_mdanalysis(u)

# Set the size of the viewer
view._remote_call('setSize', target='Widget', args=['600px', '400px'])

# Generate HTML for the viewer
with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.html') as tmp:
    nv.write_html(tmp.name, view)
    with open(tmp.name, 'r') as f:
        html_content = f.read()

# Display the viewer
st.components.v1.html(html_content, height=500, width=700)

# Frame selection slider
frame = st.slider("Select frame", 0, num_frames-1, 0)

# Update the frame
view.frame = frame

# Button to play/pause animation
if st.button("Play/Pause"):
    view._iplayer._play()