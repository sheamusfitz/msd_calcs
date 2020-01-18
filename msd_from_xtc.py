
import sys
print('Version ',sys.version)

import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda
# import pandas as pd
import re
# import json


usemda = False
if usemda == True:
    u = mda.Universe('system.gro')
    u.atoms.masses = 1

    u = mda.Merge(u.select_atoms('not (name W or name CL-)'))
    u.load_new('nojump.xtc')

    # print(u.trajectory.ts.positions[0][0])

    kalpt_1 = u.atoms[0:110]
    print(kalpt_1.positions.shape)

with open('system.top') as f:
    lines = f.readlines()
    i = -1
    while i < len(lines):
        i += 1
        if re.search(r"\[\s*molecules", lines[i]):
            print(i,lines[i])
            break
    while i < len(lines):
        i += 1
        if re.search(r"^;", lines[i]) == None:
            print(lines[i])
            print(lines[i].split())
            break
