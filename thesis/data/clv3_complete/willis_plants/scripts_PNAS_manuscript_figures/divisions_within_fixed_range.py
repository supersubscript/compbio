import cPickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as pl
from scipy.stats.stats import pearsonr

##################################################
from config_2015_02_16_plant1 import rootPath, ALTDataPathList, segDataPathList
from indices_parameters_2015_02_16_plant1 import indicesL1Init, L1centre, diceThreshold
radius = 30 
zMin = 8
################################################

region = indicesL1Init[0] #check parameters file that this corresponds to radius
indicesTStart = indicesL1Init[1]
dataLineages, scores = [], []
segData = []
cellBirths, cellDivisions, cellDivisionsPositions = [], [], [] 
#record cell cycles at pairs [(j0, m0, d0), ...], [(j1,m1,d11, d12), ...] where cell is born between time[j0] with mother m0 and time[j0+1] with daughter d0, and divides between time[j1] and time[j1+1] with mother m1 and daughters d11, d12

for i in range(len(ALTDataPathList)):
    fobj = file(ALTDataPathList[i])
    dataLineage, score = cPickle.load(fobj)
    fobj.close()
    dataLineages.append(dataLineage); scores.append(score)
    fobj = file (segDataPathList[i])    
    seg =  cPickle.load(fobj)
    fobj.close() 
    segData.append(seg)

fobj = file (segDataPathList[-1])    
seg =  cPickle.load(fobj)
fobj.close() 
segData.append(seg)

##new code from here. THIS BLOCK OF CODE WAS NOT IN ORIGINAL PNAS RESULTS. INSTEAD THE COMMENTED LINE WAS USED. THIS ENSURES THAT ALL CELLS WERE AND CHANGES RESULTS ONLY VERY MARGINALLY.
selectedIndices = []
for index in indicesTStart:
    xSq = (segData[0][0]['barycenter'][index][0] - L1centre[0][0]) * (segData[0][0]['barycenter'][index][0] - L1centre[0][0])
    ySq = (segData[0][0]['barycenter'][index][1] - L1centre[0][1]) * (segData[0][0]['barycenter'][index][1] - L1centre[0][1])
    zSq = (segData[0][0]['barycenter'][index][2] - L1centre[0][2]) * (segData[0][0]['barycenter'][index][2] - L1centre[0][2])
    if( np.sqrt(xSq + ySq + zSq) < radius and (1 in segData[0][0]['neigbourhood'][index])): 
        selectedIndices.append(index)
#########################################

#selectedIndices = indicesTStart

for i in range(0, len(segDataPathList)-1):
    allDaughters = [] 
    newBirths = []
    for mother in selectedIndices:
        daughters = [x[0] for x in dataLineages[i] if x[1] == mother]
        if len(daughters) == 2:
            newBirths.append((i, mother, daughters[0])); newBirths.append((i, mother, daughters[1])) 
        allDaughters.append(daughters)

    #for newBirths, track each cell forward to find corresponding division; if there is a division, append to allBirths, allDivisions 
    for newborn in newBirths:
         j = 1; daughters = [newborn[2]]
         while len(daughters) == 1 and i + j < len(dataLineages):
             mother = daughters[0] 
             nextDaughters = [x[0] for x in dataLineages[i + j] if x[1] == mother]; nextScore = [scores[i+j][k] for k in range(len(dataLineages[i+j])) if dataLineages[i+j][k][1] == mother]
             daughters = nextDaughters  
             j  = j + 1
         
         if(len(daughters) == 2 and min(nextScore) > diceThreshold): 
             if(1 in segData[i + 1][0]['neigbourhood'][newborn[2]] and segData[i + 1][0]['barycenter'][newborn[2]][2] > zMin and 1 in segData[i + j ][0]['neigbourhood'][daughters[0]] and 1 in segData[i + j ][0]['neigbourhood'][daughters[1]] and segData[i + j][0]['barycenter'][daughters[0]][2] >  zMin and segData[i + j][0]['barycenter'][daughters[1]][2] > zMin):
                 cellBirths.append(newborn); cellDivisions.append((i + j - 1, mother, daughters[0], daughters[1])) 
                 xSq0 = (segData[i + j][0]['barycenter'][daughters[0]][0] - L1centre[i+j][0])* (segData[i + j][0]['barycenter'][daughters[0]][0] - L1centre[i+j][0])
                 ySq0 = (segData[i + j][0]['barycenter'][daughters[0]][1] - L1centre[i+j][1])* (segData[i + j][0]['barycenter'][daughters[0]][1] - L1centre[i+j][1])
                 zSq0 = (segData[i + j][0]['barycenter'][daughters[0]][2] - L1centre[i+j][2])* (segData[i + j][0]['barycenter'][daughters[0]][2] - L1centre[i+j][2])
                 xSq1 = (segData[i + j][0]['barycenter'][daughters[1]][0] - L1centre[i+j][0])* (segData[i + j][0]['barycenter'][daughters[1]][0] - L1centre[i+j][0])
                 ySq1 = (segData[i + j][0]['barycenter'][daughters[1]][1] - L1centre[i+j][1])* (segData[i + j][0]['barycenter'][daughters[1]][1] - L1centre[i+j][1])
                 zSq1 = (segData[i + j][0]['barycenter'][daughters[1]][2] - L1centre[i+j][2])* (segData[i + j][0]['barycenter'][daughters[1]][2] - L1centre[i+j][2])
                 cellDivisionsPositions.append(0.5* (np.sqrt(xSq0 + ySq0 + zSq0) + np.sqrt(xSq1 + ySq1 + zSq1)))
    #includes all daughters
    
    #keep only indices within certain distance of centre
    fobj = file (segDataPathList[i + 1])
    dataSegT1 = cPickle.load(fobj)
    fobj.close() 
    selectedIndices = []
    for index in [item for sublist in allDaughters for item in sublist]:
        xSq = (segData[i+1][0]['barycenter'][index][0] - L1centre[i+1][0]) * (segData[i+1][0]['barycenter'][index][0] - L1centre[i+1][0])
        ySq = (segData[i+1][0]['barycenter'][index][1] - L1centre[i+1][1]) * (segData[i+1][0]['barycenter'][index][1] - L1centre[i+1][1])
        zSq = (segData[i+1][0]['barycenter'][index][2] - L1centre[i+1][2]) * (segData[i+1][0]['barycenter'][index][2] - L1centre[i+1][2])
        if( np.sqrt(xSq + ySq + zSq) < radius and (1 in segData[i+1][0]['neigbourhood'][index])): 
            selectedIndices.append(index)

    #select cells to include from daughters

np.savetxt(rootPath + '/Div_t0_new_' + str(radius) + '_' + region + '.dat', cellBirths)
np.savetxt(rootPath + '/Div_t1_new_' + str(radius) + '_' + region + '.dat', cellDivisions)

print len(cellBirths)

#identify and save all sisters
sis1Birth, sis2Birth = [], []
sis1Div, sis2Div = [], [] 
for i in range(len(cellBirths)):
    for j in range(i+1, len(cellBirths)):
        if (cellBirths[j][0] == cellBirths[i][0] and cellBirths[j][1] == cellBirths[i][1] and cellBirths[j][2] != cellBirths[i][2]): 
            [n1,n2] = [i, j] if np.random.rand() < 0.5 else [j, i]
            sis1Birth.append(cellBirths[n1]); sis2Birth.append(cellBirths[n2]); sis1Div.append(cellDivisions[n1]); sis2Div.append(cellDivisions[n2])

print len(sis1Birth)
print len(sis1Div)
np.savetxt(rootPath + '/sis1_div_t0_r_%s_z_%s_%s.dat'%(radius, zMin, region), sis1Birth)
np.savetxt(rootPath + '/sis1_div_t1_r_%s_z_%s_%s.dat'%(radius, zMin, region), sis1Div)
np.savetxt(rootPath + '/sis2_div_t0_r_%s_z_%s_%s.dat'%(radius, zMin, region), sis2Birth)
np.savetxt(rootPath + '/sis2_div_t1_r_%s_z_%s_%s.dat'%(radius, zMin, region), sis2Div)


for cell1 in cellBirths:
    count = 0
    for cell2 in cellBirths:
        if(cell1[0] ==cell2[0] and cell1[2] == cell2[2]): count = count + 1
    if(count > 1): print "duplicate!"    


fig =plt.figure() 
#ax = fig.add_subplot(1,1,1)  
#bins = [0,0.2,0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2,2.2]
plt.hist(cellDivisionsPositions, alpha=0.5, color = 'c')
#plt.xlim(0, 2.2)
plt.xlabel('position ($\mu$m)')
#plt.ylabel('frequency')
#adjustFigAspect(fig,aspect=0.8)
#fig.set_size_inches(5.5, 4.0)
#plt.gcf().subplots_adjust(bottom=0.12, left = 0.12)
#plt.yticks(np.arange(0.0, 401, 100))
plt.savefig(rootPath + "/freq_division_positions.pdf", format='pdf', dpi=300)
plt.close()
