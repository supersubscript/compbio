import cPickle
import numpy as np
import os
import sys

# with open(inFile, "rb") as input_file:
# fobj = file(sys.argv[1])
path = "/home/henrik/compbio/thesis/data/newData/plant2/nuclear_data/t0_BOA_volumes_segmentationLabels.pkl"
fobj = file(path)
data = cPickle.load(fobj)

for key, value in data.items():
  print key, value

# print data.keys()
# for mapping, score in zip(data[0], data[1]):
  # print mapping[0], mapping[1], score

# fobj.close()

### Print to file
# data = all_data["neigbourhood"]
# with open(inFile[0:(-4)] + '_neighbourhood.dat', 'a') as output_file:
#   for key, value in data.items():
#     output_file.write(str(key) + "\t" + '\t'.join(map(str, value)) + "\n")
#     # print key, ' '.join(map(str, value))
# 
# data = all_data["volumes"]
# with open(inFile[0:(-4)] + '_volumes.dat', 'a') as output_file:
#     output_file.write(str(key) + "\t" + str(value) + "\n")
# 
# data = all_data["L1"]
# with open(inFile[0:(-4)] + '_L1.dat', 'a') as output_file:
#   for cellId in data:
# #     # print cellId
#     output_file.write(str(cellId) + "\n")
#  
# data = all_data["barycenter"]
# with open(inFile[0:(-4)] + '_barycenter.dat', 'a') as output_file:
#   for key, value in data.items():
#     # print key, ' '.join(map(str, value))
#     output_file.write(str(key) + '\t' + '\t'.join(map(str, value)) + "\n")


# dictionary keys in the file:
# 'background_neighbors', 'labels', 'stem', 'decimatedSurfaceArea', 'wall_surface', 'barycenter', 'volumes', 'L1', 'neigbourhood'
