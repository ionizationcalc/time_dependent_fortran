#!/usr/bin/env python
# Plots data in unformatted outputs.
# Usage: python plot.py ../path_to_/Test_001.dat

import numpy as np
import sys
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.io import FortranFile

#
# Problem parameters
#
datafile = sys.argv[1]

#
# Read binary file
#
try:
  file = FortranFile(datafile, 'r')
except:
  print("Usage: python example-plot.py <datefile>")
  raise SystemExit

# Read the line number
nline = file.read_reals(dtype=np.int32)

# Set pdf pages
with PdfPages('./nei_example-plot.pdf') as pdf:

  for iline in range(nline[0]):
    nei_conce = file.read_reals(dtype=np.float64).reshape(30, 30)
    ei_conce = file.read_reals(dtype=np.float64).reshape(30, 30)
    
    plt.figure(figsize=(11/1.5, 8.5/1.5))
    plt.yscale("log")
    plt.ylim((1.0e-5, 1.0))
    plt.xlabel("Charge States")
    plt.ylabel("Ion Fractions")
    
    for ielem in range(30):
      natom = ielem + 1
      nstate = natom + 1
      nei_current = nei_conce[ielem]
      ei_current = ei_conce[ielem]
      if (np.sum(ei_current) != 0):
        print(f"natom={natom}, ei={ei_current[0:nstate]}")
        charge_state = np.linspace(1, nstate, nstate)
        plt.plot(charge_state, ei_current[0:nstate], ls='dotted')
        natom_str = "Natom={:d}".format(natom)
        plt.plot(charge_state, nei_current[0:nstate], label=natom_str)
    plt.legend()
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()