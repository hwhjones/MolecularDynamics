import streamlit as st
import MDAnalysis as mda
import molviewspec as mvs

# Load the multi-frame PDB file
u = mda.Universe('converted.pdb')

# Create a MolViewSpec instance
builder = mvs.create_builder()

(
    builder.load(url='converted.pdb')
    .parse(format='pdb')
    .assembly_structure(assembly_id='1')
    .component()
    .representation()
    .color(color='#1b9e77')
)
# Get the state of the builder
state = builder.get_state()

# Display the structure
st.components.v1.html(mvs.to_html(state), height=600)

# Print the state (optional, for debugging)
st.write(state)