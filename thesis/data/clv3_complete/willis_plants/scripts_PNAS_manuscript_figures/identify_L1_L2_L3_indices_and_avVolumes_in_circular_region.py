import cPickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr


############*****************###################
datafile ="/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant1/segmentation_data/0hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl"

from indices_parameters_2015_02_16_plant1 import centre3DBegin, centre3DEnd


[xRes, yRes, zRes] = (0.2635765, 0.2635765, 0.26)
###########selected circle with centre and radius given below, probably central zone.

[centreX, centreY, centreZ] = centre3DEnd
radius = 30
zThr = 10 
celldiff = 2  
# zThr. Cells that are labelled as L1 but are below this z-threshold are removed (thereby removing bottom layer)
# celldiff. L1 cells that have z-barycenter different from mean of neighbours' z-barycenter by this amount are removed 
######*************************#################

xThr = xRes * 5
yThr = yRes * 5
xMax = xRes * 512
yMax= yRes * 512

fobj = file(datafile)
dataSeg = cPickle.load(fobj)
fobj.close()

###First identify a list of cells in the L1 and L2
#Cells with barycenter below z-threshold are excluded from L1.  Also cells lying on outside of x-y frame
#This is because all cells adjacent to background are included in L1, including cells on underside.
 
L1indices = []
for index in dataSeg[0]['L1'] :   
    if((dataSeg[0]['barycenter'][index][2] > zThr) and
    (dataSeg[0]['barycenter'][index][1] > yThr) and (dataSeg[0]['barycenter'][index][1] < yMax - yThr) and
    (dataSeg[0]['barycenter'][index][0] > xThr) and (dataSeg[0]['barycenter'][index][0] < xMax - xThr)):   
        L1indices.append(index)

#remove L2 cells that are below L1 cells missing from the segmentation
incorrectIndices = []
for index in L1indices:
    array = []
    for indexNeigh in dataSeg[0]['neigbourhood'][index]:    
        if ((indexNeigh in L1indices) and indexNeigh != 1):
            array.append(dataSeg[0]['barycenter'][indexNeigh][2])    
    averageZ =  np.mean(array)
    if (dataSeg[0]['barycenter'][index][2]**2 < (averageZ - celldiff)**2):
        incorrectIndices.append(index)
print "Number of incorrect indices", len(incorrectIndices)
for index in incorrectIndices:
    L1indices.remove(index)                   

L2indices = []
for index in L1indices:
    for indexNeigh in dataSeg[0]['neigbourhood'][index]:
        if(indexNeigh != 1 and (dataSeg[0]['barycenter'][indexNeigh][1] > yThr) and (dataSeg[0]['barycenter'][indexNeigh][1] < yMax - yThr) and
        (dataSeg[0]['barycenter'][indexNeigh][0] > xThr) and (dataSeg[0]['barycenter'][indexNeigh][0] < xMax - xThr)):
            if not((1 in dataSeg[0]['neigbourhood'][indexNeigh])):   
                if (indexNeigh not in L1indices) and (indexNeigh not in L2indices): 
                    L2indices.append(indexNeigh)

L3indices = []
for index in L2indices:
    for indexNeigh in dataSeg[0]['neigbourhood'][index]:
        if(indexNeigh != 1 and (dataSeg[0]['barycenter'][indexNeigh][1] > yThr) and (dataSeg[0]['barycenter'][indexNeigh][1] < yMax - yThr) and
        (dataSeg[0]['barycenter'][indexNeigh][0] > xThr) and (dataSeg[0]['barycenter'][indexNeigh][0] < xMax - xThr)):
            if not((1 in dataSeg[0]['neigbourhood'][indexNeigh])):   
                if (indexNeigh not in L1indices) and (indexNeigh not in L2indices)  and (indexNeigh not in L3indices): 
                    L3indices.append(indexNeigh)

#for cells in the L1 and L2, include them in region of interest if they are inside the circle:
L1indicesInCircle = []
for index in L1indices:
    if((dataSeg[0]['barycenter'][index][2] - centreZ)*(dataSeg[0]['barycenter'][index][2] - centreZ) + (dataSeg[0]['barycenter'][index][1] - centreY)*(dataSeg[0]['barycenter'][index][1] - centreY) + (dataSeg[0]['barycenter'][index][0] - centreX)*(dataSeg[0]['barycenter'][index][0] - centreX) < radius*radius):
        L1indicesInCircle.append(index)

print datafile
print "L1 indices in circle : " 
print L1indicesInCircle

L2indicesInCircle = []
for index in L2indices:
    if((dataSeg[0]['barycenter'][index][2] - centreZ)*(dataSeg[0]['barycenter'][index][2] - centreZ) + (dataSeg[0]['barycenter'][index][1] - centreY)*(dataSeg[0]['barycenter'][index][1] - centreY) + (dataSeg[0]['barycenter'][index][0] - centreX)*(dataSeg[0]['barycenter'][index][0] - centreX) < radius*radius):
        L2indicesInCircle.append(index)

print "L2 indices in circle : " 
print L2indicesInCircle

L3indicesInCircle = []
for index in L3indices:
    if(((dataSeg[0]['barycenter'][index][2] - centreZ)*(dataSeg[0]['barycenter'][index][2] - centreZ) + dataSeg[0]['barycenter'][index][1] - centreY)*(dataSeg[0]['barycenter'][index][1] - centreY) + (dataSeg[0]['barycenter'][index][0] - centreX)*(dataSeg[0]['barycenter'][index][0] - centreX) < radius*radius):
        L3indicesInCircle.append(index)

print "L3 indices in circle : " 
print L3indicesInCircle

volumes = []
for index in (L1indicesInCircle + L2indicesInCircle):      
        volumes.append(dataSeg[0]['volumes'][index])
print np.mean(volumes)

#output = open(folder + datafile[:-4] + "indicesCentralZone.txt", "w")
#output.write("Cell indices of L1 in circular region, (centreX, centreY)  = (" + str(centreX) + "," + str(centreY) +  "), radius = " + str(radius) + "\n")
#output.write(str(L1indicesInCircle) + "\n")

#output.write("Cell indices of L2 in circular region, (centreX, centreY)  = (" + str(centreX) + "," + str(centreY) +  "), radius = " + str(radius) + "\n")
#output.write(str(L2indicesInCircle) + "\n")  

#output.write("Cell indices of L3 in circular region, (centreX, centreY)  = (" + str(centreX) + "," + str(centreY) +  "), radius = " + str(radius) + "\n")
#output.write(str(L3indicesInCircle) + "\n")              
  

#output.write("Cell indices of L1 and L2 in circular region, (centreX, centreY)  = (" + str(centreX) + "," + str(centreY) +  "), radius = " + str(radius) + "\n")
#output.write(str(L1indicesInCircle + L2indicesInCircle) + "\n")

#output.write("Mean volume in region: \n")
#output.write(str(np.mean(volumes)))
#output.close()
