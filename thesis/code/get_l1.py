# get.l1
import cPickle
import os

plants = [1, 2, 4, 13, 15, 18]
data_path = "../data/newData/"


# Get absolute paths to directories
def absoluteFilePaths(directory):
   for dirpath,_,filenames in os.walk(directory):
       for f in filenames:
           yield os.path.abspath(os.path.join(dirpath, f))


for plant in plants:
  files = absoluteFilePaths(data_path + "plant" + str(plant) + "/segmentation_data/")
  files = filter(lambda x: x.endswith('.pkl'), files)
  
  for inFile in files:
    # Read in data
    with open(inFile, "rb") as input_file:
      fobj = file(inFile)
      all_data = cPickle.load(fobj)
      all_data = all_data[0]
      
      ### Print to file
      # data = all_data["neigbourhood"]
      # with open(inFile[0:(-4)] + '_neighbourhood.dat', 'a') as output_file:
      #   for key, value in data.items():
      #     output_file.write(str(key) + "\t" + '\t'.join(map(str, value)) + "\n")
      #     # print key, ' '.join(map(str, value))
      # 
      # data = all_data["volumes"]
      # with open(inFile[0:(-4)] + '_volumes.dat', 'a') as output_file:
      #   for key, value in data.items():
      # #     # print key, value
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
      
      fobj.close()

# dictionary keys in the file:
# 'background_neighbors', 'labels', 'stem', 'decimatedSurfaceArea', 'wall_surface', 'barycenter', 'volumes', 'L1', 'neigbourhood'
