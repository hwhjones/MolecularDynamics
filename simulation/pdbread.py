from ase.io import read
from ase.io import write as ase_write

# Read Traj file
traj = read('nacl_water.traj', index=":")
ase_write('nacl_water.pdb', traj, format='proteindatabank')