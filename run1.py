import streamlit as st
import nglview as nv
import nglview
from ase.io import read
import io
import numpy as np
from PIL import Image

# Read the trajectory file
traj = read('nacl_water.traj', index=':')

# Create a Streamlit app
st.title("Molecular Dynamics Trajectory Viewer")

# Create an NGLView widget with the trajectory
view = nv.show_asetraj(traj)

nglview.write_html("index.html", [view])

# Read the contents of the file
with open("index.html", 'r') as file:
    html_content = file.read()

st.components.v1.html(html_content, height=620, width=620)