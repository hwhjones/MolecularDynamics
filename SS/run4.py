from Bio import PDB
from io import StringIO
import streamlit as st
import py3Dmol
from stmol import showmol
import nglview as nv

def read_pdb(file_path):
    with open(file_path, 'r') as file:
        return file.read()

def render_mol(pdb_string, frame=0):
    # Create a PDB parser
    parser = PDB.PDBParser(QUIET=True)
    
    # Parse the PDB string
    structure = parser.get_structure("molecule", StringIO(pdb_string))
    
    # Get the specified model (frame)
    model = structure[frame + 1]  # PDB model numbering starts at 1
    
    # Convert the model to PDB string
    io = PDB.PDBIO()
    io.set_structure(model)
    frame_pdb = StringIO()
    io.save(frame_pdb)
    frame_pdb_string = frame_pdb.getvalue()
    
    # Create and setup the view
    view = py3Dmol.view(width=800, height=400)
    view.addModel(frame_pdb_string, 'pdb')
    view.setStyle({'sphere': {'color': 'orange', 'radius': 0.7}})
    view.zoomTo()
    view.setBackgroundColor('white')
    
    # Show the molecule
    showmol(view, height=400, width=800)

# Streamlit app
st.title('Copper Atoms Trajectory Viewer')

# Read the PDB file
pdb_file_path = 'converted.pdb'  # Replace with your actual file path
pdb_string = read_pdb(pdb_file_path)

# Count the number of frames
structure = PDB.PDBParser(QUIET=True).get_structure("molecule", StringIO(pdb_string))
num_frames = len(structure)

# Frame selection
frame = st.slider('Select frame', 0, num_frames - 1, 0)

# Render the molecule
st.subheader("Molecule Visualization")
mol_html = render_mol(pdb_string, frame)
st.components.v1.html(mol_html, height=400, width=800)

# Display frame information
st.write(f"Showing frame {frame + 1} of {num_frames}")

# Display atom coordinates for the current frame
st.subheader("Atom Coordinates")
model = structure[frame + 1]
for atom in model.get_atoms():
    st.text(f"{atom.get_full_id()[3][1]:>5} {atom.get_name():>4} {atom.get_coord()[0]:>8.3f} {atom.get_coord()[1]:>8.3f} {atom.get_coord()[2]:>8.3f}")

# Display raw PDB content
st.subheader("Raw PDB Content")
st.text(pdb_string)