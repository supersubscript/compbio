#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May 21 18:36:36 2017

@author: henrik
"""
import os
import pandas as pd
import numpy as np
os.chdir("/home/henrik/compbio/thesis/code/")

file = open("../data/rk.test", "r")
data = pd.read_csv("../data/rk.test", sep="|", names=['A'], index_col = False)
cols = list(range(0,18))
data[cols] = data.A.str.split('\t', n=18, expand=True)
data = data.drop("A", 1)
data = data.astype(np.float)

# Extract the fun stuff
nDatapoints = int(data[0][0])
nCells = int(data[0][1])
nSpecies = int(data[1][1])

data = data[np.isfinite(data[3])] # Take out rows with NaNs
data.index = [ ii for ii in xrange(0,len(data)) ] # Reset indices

trajectories = np.zeros((int(nSpecies), int(nCells), int(nDatapoints)), dtype = np.float)

for ii in xrange(0, nDatapoints):
    for jj in xrange(0, nSpecies):
        trajectories[jj, :, ii] = data[jj][ii * nCells:(ii + 1) * nCells]

##########
import numpy as np
import matplotlib.pyplot as plt


# Example data
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
for ii in xrange(0, nCells):
  plt.plot(trajectories[6,ii,:100])

plt.xlabel(r'Time',fontsize=16)
plt.ylabel(r'Expression',fontsize=16)
#plt.title(r"\TeX\ is Number "
#          r"$\displaystyle\sum_{n=1}^\infty\frac{-e^{i\pi}}{2^n}$!",
#          fontsize=16, color='gray')

# Make room for the ridiculously large title.
plt.subplots_adjust(top=0.8)

#plt.savefig('tex_demo')
plt.show()
