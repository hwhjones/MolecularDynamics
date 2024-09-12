import torch
from ase.io import read
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io.trajectory import Trajectory
from ase.io import write as ase_write
from ase.md import MDLogger
from ase import units
from orb_models.forcefield import atomic_system
from orb_models.forcefield import pretrained
from ase.calculators.calculator import Calculator, all_properties
import numpy as np
import os
from ase.md.langevin import Langevin  # Langevin thermostat for NVT

# Define input file name and cell size at the top
input_file = "nacl_water.xyz"
cell_size = 13.92  # in Angstroms

# Generate output file name (change extension to .pdb)
output_file = os.path.splitext(input_file)[0] + ".pdb"

# Define a custom Orb-d3-v1 calculator with CUDA support
class OrbD3Calculator(Calculator):
    implemented_properties = ['energy', 'forces', 'stress']

    def __init__(self, model, **kwargs):
        super().__init__(**kwargs)
        self.model = model

        # Check if CUDA is available and move the model to GPU if possible
        if torch.cuda.is_available():
            self.model.cuda()
        else:
            print("CUDA is not available, running on CPU")

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_properties):
        # Ensure the atoms object is set correctly
        super().calculate(atoms, properties, system_changes)

        # Convert the ASE atoms to a graph representation
        graph = atomic_system.ase_atoms_to_atom_graphs(atoms)

        # Move the graph to GPU if CUDA is available
        if torch.cuda.is_available():
            graph = graph.to('cuda')  # Use .to('cuda') to move the graph to GPU

        # Perform the prediction using orb-d3-v1
        result = self.model.predict(graph)

        # Store results in self.results dictionary for ASE compatibility
        self.results['energy'] = result["graph_pred"].detach().cpu().numpy().sum()  # Sum total energy
        self.results['forces'] = result["node_pred"].detach().cpu().numpy()  # Forces
        self.results['stress'] = result["stress_pred"].detach().cpu().numpy()  # Stress (if needed)

        # After prediction, print the forces
        #print("Predicted forces (eV/Å):", self.results['forces'])

# Load the checkpoint manually
#local_model_path = "/scratch/tn51/ttd110/orb/orb-d3-v1-20240902.ckpt"

# Load the model using the local checkpoint path
#orbff = pretrained.orb_d3_v1(weights_path=local_model_path)
orbff = pretrained.orb_d3_v1()

# Load XYZ file into ASE Atoms object
atoms = read(input_file)

# Update cell size
atoms.set_cell([cell_size, cell_size, cell_size])
atoms.set_pbc([True, True, True])  # Apply periodic boundary conditions

# Attach the custom Orb-d3-v1 calculator to the atoms object
atoms.calc = OrbD3Calculator(model=orbff)

# Set up the initial velocities corresponding to a temperature of 300 K
MaxwellBoltzmannDistribution(atoms, temperature_K=300)

# Define the Langevin thermostat for NVT
temperature_K = 300  # Temperature in Kelvin
friction = 0.05  # Friction coefficient, adjust as needed
timestep = 0.5 * units.fs

dyn = Langevin(atoms, timestep, temperature_K=temperature_K, friction=friction)

# Function to write PDB file
def write_pdb(atoms=atoms):
    ase_write(output_file, atoms, format='pdb', append=True)

# Attach the PDB writing function instead of Trajectory
dyn.attach(write_pdb, interval=20)  # Save every 20 steps

# Remove or comment out the Trajectory-related lines
# trajectory = Trajectory(output_file, "w", atoms)
# dyn.attach(trajectory.write, interval=20)

# Optional: Add an MDLogger to print energy and forces
dyn.attach(MDLogger(dyn, atoms, "md_nvt.log", header=True, stress=False, peratom=True), interval=10)

# Run the simulation for 1,000,000 time steps
dyn.run(1000000)