from ase.io import read

traj = read('nacl_water.traj', index=':')

from ase.visualize import view

view(traj)