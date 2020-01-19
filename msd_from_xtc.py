
import sys
print('Version ',sys.version)

import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda
import re
import json
import readline


def indexer():
  """
  This function (at the current moment) runs on the current folder and creates an index matrix of the names and numbers of molecules
  """
  molecules = []
  with open('system.top') as f:
    lines = f.readlines()
    in_molecules = False
    for line in lines:
      if in_molecules == True:
        if not re.search(r"^;", line):
          # print(line)
          molecules.append(line.split())
          # print(molecules)
      if re.search(r"\[\s*molecules", line):
        # print(line)
        in_molecules = True
  for molecule in molecules:
    molecule[1] = int(molecule[1])
  # print(molecules)
  # at this point, the 'molecules' object lists the names (conflicts possible) of the molecules and the quantity of them
  molecules_of_interest = []
  for molecule in molecules:
    interest_input = input('Should I analyze \033[1m'+molecule[0]+'\033[0m? (y/n/q) :')
    if interest_input == 'y':
      molecules_of_interest.append(molecule)
    elif interest_input == 'q':
      break
  # print(molecules_of_interest)

  with open("mols_of_interest.json", "w") as f:
    json.dump(dict(molecules_of_interest), f)

def msd():
  with open('./mols_of_interest.json') as f:
    data = json.load(f)
  u = mda.Universe('system.gro')
  u.atoms.masses = 1
  u = mda.Merge(u.select_atoms('not (name W or name CL-)'))
  u.load_new('nojump.xtc')
  print(u.trajectory.coordinate_array("afc"))
  print(sum(data.values()))
  for species in data:
    # print(species, data[species])
    print()

    for molecule in range(data[species]):
      # print(molecule)
      print()









if __name__ == "__main__":
  try:
    with open('./mols_of_interest.json') as f:
      print('yeah it\'s there')
  except IOError:
    indexer()
  msd()
