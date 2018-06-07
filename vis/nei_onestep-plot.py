#!/usr/bin/env python
# Plots data in unformatted outputs.
# Usage: python nei_onestep-plot.py ../path_to_/test_onestep_ei.dat

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
  f = FortranFile(datafile, 'r')
except:
  print("Usage: python nei_onestep-plot.py <datefile>")
  raise SystemExit
# Read Te and ne
[te_sta, te_end, ne] = f.read_reals(dtype=np.float64)

# Read the ion fractions at the beginning and end time
fraction_sta_ei = f.read_reals(dtype=np.float64).reshape(30, 30)
fraction_end_ei = f.read_reals(dtype=np.float64).reshape(30, 30)
print(fraction_sta_ei[1, :])
print(fraction_sta_ei[25, :])
print()

# Read number of cases
ntime = f.read_ints()
print(f"ntime={ntime}")

# Set pdf pages
with PdfPages('./nei_onestep-plot.pdf') as pdf:
  plt.figure(figsize=(11/1.5, 8.5/1.5))
  plt.yscale("log")
  plt.ylim((1.0e-5, 1.0))
  plt.xlabel("Charge States")
  plt.ylabel("Iron Fractions")
  plt.title('Te={:10.2e}K, ne={:10.2e}cm^-3'.format(te_end, ne))
  plt.grid(True)

  # Plot Fe_EI at the beginning and end point
  natom = 26
  nstate = natom + 1
  charge_state = np.linspace(1, 27, 27)
  plt.plot(charge_state, fraction_sta_ei[natom-1, 0:nstate],c='b', 
    label='start')
  plt.plot(charge_state, fraction_sta_ei[natom-1, 0:nstate],'o',c='b')
  plt.plot(charge_state, fraction_end_ei[natom-1, 0:nstate], c='r', 
    label='EI end (Te={:7.2e}K)'.format(te_end))
  plt.plot(charge_state, fraction_end_ei[natom-1, 0:nstate], 'o', c='r')

  for itime in range(ntime[0]):
    dt = f.read_reals(dtype=np.float64)
    print(f"dt={dt}")
    fraction_nei = f.read_reals().reshape(30, 30)
    label_str = 'NEI t={:.1f}s'.format(dt[0])
    plt.plot(charge_state, fraction_nei[natom-1, 0:nstate], ls='--', 
            label=label_str)
  plt.legend()
  pdf.savefig()  # saves the current figure into a pdf page
  plt.close()