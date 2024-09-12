from ase.io import read

traj = read('nacl_water.traj', index=':')

from ase.visualize import view

view(traj)

import pandas as pd
import numpy as np

df = pd.read_csv('md_nvt.log')

print (df.head)