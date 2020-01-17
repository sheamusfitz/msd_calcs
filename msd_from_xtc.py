import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda
# import pandas as pd
# import re
# import json

u = mda.Universe("system.gro")
u.atoms.masses = 1

u2 = mda.Merge(u.select_atoms("not (name W or name CL-)"))
u2.load_new("nojump.xtc")

# print(u.trajectory.ts.positions[0][0])

kalpt_1 = u2.atoms[0:110]
# TODO I need to read in this "110" from main file
# I have to loop over the thing to get this to a numpy array


print(kalpt_1.positions.shape)
