import streamlit as st
import py3Dmol
from stmol import showmol

def read_pdb(file_path):
    with open(file_path, 'r') as file:
        return file.read()

def render_mol(pdb_string, frame=0):
    pdb_lines = pdb_string.split('\n')
    start_index = pdb_lines.index(f'MODEL     {frame + 1}') + 1
    end_index = pdb_lines.index('ENDMDL', start_index)
    frame_pdb = '\n'.join(pdb_lines[start_index:end_index])
    print(frame_pdb)
    print(start_index)
    print(end_index)
    view = py3Dmol.view(width=800, height=400)
    view.addModel(frame_pdb, 'pdb', {'keepH': True})
    view.setBackgroundColor('white')

    # Set style for all atoms as spheres
    view.setStyle({}, {'sphere': {'radius': 0.3}})

    # Adjust sphere size and color for different elements
    view.addStyle({'atom':'O'}, {'sphere': {'radius': 0.66, 'color': 'red'}})
    view.addStyle({'atom':'H'}, {'sphere': {'radius': 0.31, 'color': 'black', 'opacity':0.8}})
    view.addStyle({'atom':'Na'}, {'sphere': {'radius': 1.02, 'color': 'purple'}})
    view.addStyle({'atom':'Cl'}, {'sphere': {'radius': 0.99, 'color': 'green'}})

    view.zoomTo()

    showmol(view, height=400, width=800)

# Streamlit app
st.title('NaCl Atoms Trajectory Viewer')

# Read the PDB file
pdb_file_path = 'nacl_water.pdb'  # Replace with your actual file path
pdb_string = read_pdb(pdb_file_path)

# Count the number of frames
num_frames = pdb_string.count('MODEL')
st.write(f"Number of frames detected: {num_frames}")

# Frame selection
frame = st.slider('Select frame', 0, num_frames - 1, 0)

 # Render the molecule
st.subheader("Molecule Visualization")
render_mol(pdb_string, frame)

# Display frame information
st.write(f"Showing frame {frame + 1} of {num_frames}")

# Display atom coordinates for the current frame
st.subheader("Atom Coordinates")
frame_lines = pdb_string.split('MODEL')[frame + 1].split('\n')
for line in frame_lines:
    if line.startswith('ATOM'):
        st.text(line)

# Display raw PDB content
st.subheader("Raw PDB Content")
st.text(pdb_string)