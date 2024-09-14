import streamlit as st
import py3Dmol
import time
from stmol import showmol
import base64
from io import BytesIO
from PIL import Image

pdb_file_path = 'data/nacl_water.pdb'  # Replace with your actual file path

def read_pdb(file_path):
    with open(file_path, 'r') as file:
        return file.read()

def parse_pdb(pdb_string, frame=0): # Parses PDB file into separate frames
    pdb_lines = pdb_string.split('\n')
    start_index = pdb_lines.index(f'MODEL     {frame + 1}') + 1
    end_index = pdb_lines.index('ENDMDL', start_index)
    frame_pdb = '\n'.join(pdb_lines[start_index:end_index])
    return frame_pdb

def render_mol(frame_pdb):
    view = py3Dmol.view(width=800, height=400)
    view.addModel(frame_pdb, 'pdb', {'keepH': True})
    view.setBackgroundColor('white')

    # Set style for all atoms as spheres
    view.setStyle({}, {'sphere': {'radius': 0.3}})

    # Adjust sphere size and color for different elements
    view.addStyle({'atom':'O'}, {'sphere': {'radius': 0.66, 'color': 'red'}})
    view.addStyle({'atom':'H'}, {'sphere': {'radius': 0.31, 'color': 'white'}})
    view.addStyle({'atom':'Na'}, {'sphere': {'radius': 1.02, 'color': 'purple'}})
    view.addStyle({'atom':'Cl'}, {'sphere': {'radius': 0.99, 'color': 'green'}})

    view.zoomTo()

    showmol(view, height=400, width=800)

    return view

def clear_mol(view):
    view.clear()

# Streamlit app
st.title('NaCl Atoms Trajectory Viewer')

# Read the PDB file
pdb_string = read_pdb(pdb_file_path)

# Count the number of frames
num_frames = pdb_string.count('MODEL')
# st.write(f"Number of frames detected: {num_frames}")

select = st.checkbox("Automatic Render")

if select:
    # Render the molecule
    st.subheader("Molecule Visualization")

    # Create a placeholder for the molecule viewer
    mol_viewer = st.empty()

    for frame in range(1, num_frames):
        frame_pdb = parse_pdb(pdb_string, frame)
        # Display frame information

        # Use the placeholder to display the molecule
        with mol_viewer.container():
            st.write(f"Showing frame {frame} of {num_frames}")
            view = render_mol(frame_pdb)
            # Capture the PNG image
            time.sleep(0.10)
    # Generate APNG
    apng_data = view.apng(num_frames)

    print(apng_data)


else:
    # Frame selection
    frame = st.slider('Select frame', 0, num_frames - 1, 0)

    # Render the molecule
    st.subheader("Molecule Visualization")

    frame_pdb = parse_pdb(pdb_string, frame)
    view = render_mol(frame_pdb)

    # Display frame information
    st.write(f"Showing frame {frame + 1} of {num_frames}")

    # Display atom coordinates for the current frame
    st.subheader("Atom Coordinates")
    frame_lines = frame_pdb.split('\n')
    for line in frame_lines:
        if line.startswith('ATOM'):
            st.text(line)

    # Display raw PDB content
    st.subheader("Raw PDB Content")
    st.text(pdb_string)