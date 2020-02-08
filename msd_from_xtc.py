#!python3
# import sys
# print('Version ',sys.version)

import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda
import re
import json
import readline
from termcolor import cprint
import os, sys, getopt

def dry_gro(input_file = 'system.gro'):
  with open(input_file, 'r') as input:
    with open('dry.gro', 'w') as output:
      for line in input:
        if not re.search(r"[0-9\s]{5}W", line):
          output.write(line)

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

def msd(traj_name='nojump.xtc'):
  global length2, tracking_indices, dt, u, mol_count_dict
  with open('./mols_of_interest.json') as f:
    mol_count_dict = json.load(f)
  u = mda.Universe('system.gro')
  u.atoms.masses = 1
  u = mda.Merge(u.select_atoms('not (name W or name CL-)'))
  u.load_new(traj_name)

  with open('/Users/shea/dcuments/research/msd_calcs/tracking_indices.json') as f:
    tracking_indices = json.load(f)
  # print(tracking_indices['kalp21t'])
  # print(sum(mol_count_dict.values()))
  dt = u.trajectory[1].dt *10**-12
  # print(dt)
  length2 = 0
  for species in mol_count_dict:
    print(len(mol_count_dict), mol_count_dict[species], len(u.trajectory))
    length2 = max([mol_count_dict[species], length2])
  # print(length2)

  tracked_coords = coord_array()
  msd_array = msd_calc(tracked_coords)

def coord_array():
  first_atom = 0
  tracked_coords = np.zeros((len(mol_count_dict),length2,len(u.trajectory),3), dtype=float)
  # this array is of the form:
  #   [species] x [molecule number] x [timestep] x [x,y,θ]

  i_s = 0

  for species in mol_count_dict:
    i_m = 0
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

      # plt.plot(thetas*180/np.pi)
      # plt.title("Angle of rotation (degrees)")
      # plt.show()

      tracked_coords[i_s, i_m, :, 0:2] = position_data[:, 0, 0:2]
      tracked_coords[i_s, i_m, :, 2] = thetas


      # plt.plot(tracked_coords[i_s, i_m, :, 0], tracked_coords[i_s, i_m, :, 1])

      first_atom += tracking_indices[species]['length']
      i_m += 1

    i_s += 1
  # plt.title("Center of Mass Position")
  # plt.show()

  print(tracked_coords[0,0,:,0])
  tracked_coords /= 10
  print(tracked_coords[0,0,:,0])

  # plt.plot(tracked_coords[0,0, :, 0], tracked_coords[0,0, :, 1])
  # plt.title("Center of Mass Position")
  # plt.show()

  return(tracked_coords)

def msd_calc(tracked_coords):
  msd_array = np.zeros((len(mol_count_dict), length2, int(np.log2(len(u.trajectory)))+1, 2), dtype=float)
  # the dimensions of this array are:
  #   [species] x [molecule num.] x [spacing power] x [xy, θ]
  # print((tracked_coords[0,0,0, 0:2]-tracked_coords[0,0,4, 0:2]))
  # print((tracked_coords[0,0,0, 0:2]-tracked_coords[0,0,4, 0:2])**2)
  # print(sum((tracked_coords[0,0,0, 0:2]-tracked_coords[0,0,4, 0:2])**2))

  mmsd_array = np.zeros((len(mol_count_dict), int(np.log2(len(u.trajectory)))+1, 2), dtype=float)
  for s, species in enumerate(mol_count_dict):
    for m in range(mol_count_dict[species]):
      for i in range(len(msd_array[0,0])):
        ii = 2**i

        msd_array[s, m, i, 0] = sum([sum((tracked_coords[s, m, ts, 0:2]-tracked_coords[s, m, ts+ii, 0:2])**2) for ts in range(len(u.trajectory)-ii)])/(len(u.trajectory)-ii)
        # print(msd_array[s,m,i,0])
      # plt.plot([1,2,4,8,16,32], msd_array[s,m,:,0])
    mmsd_array[s,:,0] = np.mean(msd_array[s,0:mol_count_dict[species],:,0], axis=0)
    t_arr = np.array([2**n for n in range(len(mmsd_array[s,:,0]))])*dt
    # plt.plot(t_arr*(10**6), mmsd_array[s,:,0])
  # plt.xscale('log')
  # plt.yscale('log')
  # plt.show()

  str1 = os.path.basename(os.getcwd())
  outname = 'msd'+str1
  # print(outname)
  np.save(outname, mmsd_array)

  return(mmsd_array)

if __name__ == "__main__":
  try:
    with open('./mols_of_interest.json') as f:
      print('mols_of_interest.json exists, good')
  except IOError:
    indexer()

  try:
    opts, args = getopt.getopt(sys.argv[1:],"t:")
  except getopt.GetoptError:
    print('usage: you did an option wrong i guess')
    sys.exit(2)
  dry_gro()

  if len(opts) == 0:
    filename = 'nojump.xtc'
  else:
    for opt, arg in opts:
      if opt == "-t":
        filename = arg
      else:
        filename = 'nojump.xtc'
  msd(filename)
