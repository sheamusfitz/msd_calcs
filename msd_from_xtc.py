
# import sys
# print('Version ',sys.version)

import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda
import re
import json
import readline
from termcolor import cprint

# TODO: make a function that creates a "dry" .gro file

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
    mol_count_dict = json.load(f)
  u = mda.Universe('system.gro')
  u.atoms.masses = 1
  u = mda.Merge(u.select_atoms('not (name W or name CL-)'))
  u.load_new('nojump.xtc')

  with open('./tracking_indices.json') as f:
    tracking_indices = json.load(f)
  print(tracking_indices['kalp21t'])

  print(sum(mol_count_dict.values()))

  first_atom = 0
  # print(len(u.trajectory), 'frames')
  # coms = np.zeros(
  #   (sum(mol_count_dict[s] for s in mol_count_dict), len(u.trajectory), 3)
  # )
  # rotators = np.zeros(
  #   (sum(mol_count_dict[s] for s in mol_count_dict), len(u.trajectory), 3)
  # )
  for species in mol_count_dict:
    for molecule in range(mol_count_dict[species]):
      cprint(species+'\t'+str(molecule), attrs=['bold'])

      mol = u.atoms[np.array(tracking_indices[species]['com index'])+first_atom]
      position_data = np.array([
        (
          mol.center_of_geometry(),
          u.atoms[tracking_indices[species]['rotator index']+first_atom].position
        ) for ts in u.trajectory], dtype=object)
      position_data[:,1] -= position_data[:,0]
      # print(position_data,'\n\n')
      # print(position_data[:], '\n\n')
      # print(position_data[:,0,0], '\n\n')

      thetas = np.array(
        [
          (
            np.arctan2(ts[1,1], ts[1,0])
          ) for ts in position_data
        ]
      )
      # TODO: this doesn't deal with "jumping over" from -3.14 radians to +3.14 radians
      # I'm thinking of doing some sort of loop that checks each step and goes "is theta_i+1 - theta_i > pi? -> add pi" and vice versa for the difference being less than -pi (or whatever.) I probably have signs wrong in there because math is hard

      plt.plot(thetas*180/np.pi)
      plt.title("Angle of rotation (degrees)")
      plt.show()

      first_atom += tracking_indices[species]['length']

  # This is a plot of the motion of the center of mass of the object in question
  # plt.plot(position_data[:,0,0], position_data[:,0,1])
  # plt.title("Center of Mass Position")
  # plt.show()







if __name__ == "__main__":
  try:
    with open('./mols_of_interest.json') as f:
      print('yeah it\'s there')
  except IOError:
    indexer()
  msd()
