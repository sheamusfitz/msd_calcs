# now I need to get the length of each molecule
for molecule in molecules_of_interest:
  with open(molecule+'.itp') as f:
    lines = f.readlines()
    # find [ atoms ]
    i = 0
    while i < len(lines):
      if re.search(r"\[\s*atoms", line):
        i += 1
        break
      i += 1
    while i < len(lines):
      if re.search(r"\[\s*",line):
        i -= 1
        break
      i += 1
    while i < len(lines) and i > 0:
      if re.search(r"[0-9]"):
