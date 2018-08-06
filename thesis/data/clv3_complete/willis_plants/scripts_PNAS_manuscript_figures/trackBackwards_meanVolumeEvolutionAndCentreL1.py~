import cPickle
import numpy as np
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
#############################
from config_2015_02_16_plant1 import savePath, timesList, ALTDataPathList, segDataPathList
from indices_parameters_2015_02_16_plant1 import indicesTEndBack, centre3DEnd
##############################
label = 'p1'

indicesTEnd = indicesTEndBack
region = indicesTEnd[0]
stretchingFactor = [timesList[i][3] for i in range(0,len(timesList))]
times = [float(timesList[i][0][0:-3]) for i in range(0, len(timesList))]
trackedIndices = [indicesTEnd[1]]
trackedVolumes, volumesNoErrors = [], []
meanVolumes, stdVolumes = [], []
cellClosestCentre = []
centrePositions = [] 

for i in range(len(segDataPathList)-1, 0, -1):

   #t1 file
    fobj = file (segDataPathList[i])
    dataSegT1 = cPickle.load(fobj)
    fobj.close() 

    #t0 file
    fobj = file(segDataPathList[i - 1])
    dataSegT0 = cPickle.load(fobj)
    fobj.close()

    #t0 to t1 ALT file
    fobj = file(ALTDataPathList[i- 1])
    dataLineage, scores = cPickle.load(fobj)
    fobj.close()

    print "Files :"
    print segDataPathList[i]
    print segDataPathList[i-1]
    print ALTDataPathList[i-1], "\n"

    if i == len(segDataPathList)-1: 
        print dataSegT1[0]['volumes']
        for index in indicesTEnd[1]:
            print index
            print dataSegT1[0]['volumes'][index]
        trackedVolumes.append([dataSegT1[0]['volumes'][index] for index in indicesTEnd[1]])
        meanVolumes.append( np.mean(trackedVolumes[-1]))
        stdVolumes.append( np.std(trackedVolumes[-1]))
        print "Mean volume at time ", times[i],  ": ", meanVolumes[-1]

        distances = [] 
        for index in indicesTEnd[1]:
            xSq = (dataSegT1[0]['barycenter'][index][0] - centre3DEnd[0]) * (dataSegT1[0]['barycenter'][index][0] - centre3DEnd[0])
            ySq = (dataSegT1[0]['barycenter'][index][1] - centre3DEnd[1]) * (dataSegT1[0]['barycenter'][index][1] - centre3DEnd[1])
            zSq = (dataSegT1[0]['barycenter'][index][2] - centre3DEnd[2]) * (dataSegT1[0]['barycenter'][index][2] - centre3DEnd[2])
            distances.append(np.sqrt(xSq + ySq + zSq))
        #create array of the index in indicesTStart of cells closest to the centre.
        spDist = sorted(distances)[3]
        for j in range(0, len(indicesTEnd[1])):
            if distances[j] < spDist: cellClosestCentre.append(j)

        xCentre = []
        yCentre = []
        zCentre = []
        for cellInd in cellClosestCentre:
            xCentre.append(dataSegT1[0]['barycenter'][indicesTEnd[1][cellInd]][0])
            yCentre.append(dataSegT1[0]['barycenter'][indicesTEnd[1][cellInd]][1])
            zCentre.append(dataSegT1[0]['barycenter'][indicesTEnd[1][cellInd]][2])
        centrePositions.append((np.mean(xCentre), np.mean(yCentre), np.mean(zCentre)))        

    trackedIndices.append([])
    trackedVolumes.append([])  
    volumesNoErrors.append([]) 

    ################compute 

    trackedCells = [] 
 
    for index in trackedIndices[-2]: trackedCells.append([x for x in dataLineage if x[0] == index])
    multipleParents, zeroParents, diceAv = 0,0,0 
    sumVol = 0.0 
    numCells = 0.0
    #print "Tracked cells for these files : ", trackedCells , "\n"
    for pair in trackedCells:
        if len(pair) > 1: 
            multipleParents = multipleParents + 1
            trackedIndices[-1].append(-1) 
            trackedVolumes[-1].append(-1) 
        elif len(pair) < 1: 
            zeroParents = zeroParents + 1
            trackedIndices[-1].append(-1) 
            trackedVolumes[-1].append(-1) 
        else:  
            trackedIndices[-1].append(pair[0][1])  
            trackedVolumes[-1].append(dataSegT0[0]['volumes'][pair[0][1]])
            volumesNoErrors[-1].append(dataSegT0[0]['volumes'][pair[0][1]])

    xCentre = []
    yCentre = []
    zCentre = []
    for index in cellClosestCentre: 
        xCentre.append(dataSegT0[0]['barycenter'][trackedIndices[-1][index]][0])
        yCentre.append(dataSegT0[0]['barycenter'][trackedIndices[-1][index]][1])
        zCentre.append(dataSegT0[0]['barycenter'][trackedIndices[-1][index]][2])
    centrePositions.append((np.mean(xCentre), np.mean(yCentre), np.mean(zCentre)))        
       
            
    print "Num. of cells,", "Num. with no parent,", "Num. with multiple parents : "  
    print len(trackedIndices[-1]), zeroParents, multipleParents,'\n'
    meanVolumes.append(np.mean(volumesNoErrors[-1]))
    stdVolumes.append(np.std(volumesNoErrors[-1]))

meanVolumes = meanVolumes[::-1]
stdVolumes = np.array(stdVolumes[::-1])
centrePositions = centrePositions[::-1]
noStretchcorrectedMeanVolumes = np.array([meanVolumes[i]*stretchingFactor[i] for i in range(0, len(meanVolumes))])
noStretchcorrectedStdVolumes = np.array([stdVolumes[i]*stretchingFactor[i] for i in range(0, len(stdVolumes))])
meancorrectedMeanVolumes = np.array([meanVolumes[i]/meanVolumes[i] for i in range(0, len(meanVolumes))])
meancorrectedStdVolumes = np.array([stdVolumes[i]/meanVolumes[i] for i in range(0, len(stdVolumes))])


print times, len(times)
print list(set(trackedIndices[-1]))
print "L1centre = ", centrePositions
print "meanVolumes = ", meanVolumes, len(meanVolumes)
print "mean volumes, no stretching factor = ", noStretchcorrectedMeanVolumes, len(noStretchcorrectedMeanVolumes)

meanVolumes = np.array(meanVolumes)

plt.plot(np.array(times), meanVolumes, 'b', np.array(times), meanVolumes)
plt.fill_between(np.array(times), meanVolumes - stdVolumes, meanVolumes + stdVolumes, facecolor = 'blue', alpha = 0.1)
plt.xlabel('time (hrs)')
plt.ylabel('cell mean volume with stretching correction')
plt.title(region)
plt.savefig(savePath + "%s_cell_mean_volume_stretch_correction_%s.pdf"%(label, region), format='pdf', dpi=300)
plt.close()

plt.plot(np.array(times), noStretchcorrectedMeanVolumes, 'ro', np.array(times), noStretchcorrectedMeanVolumes)
plt.plot(np.array(times), noStretchcorrectedMeanVolumes + noStretchcorrectedStdVolumes, 'r--')
plt.plot(np.array(times), noStretchcorrectedMeanVolumes - noStretchcorrectedStdVolumes, 'r--')
plt.xlabel('time (hrs)')
plt.ylabel('cell mean +/- std volume no stretch correction')
plt.title(region)
plt.savefig(savePath + "%s_cell_mean_volume_no_stretch_correction_%s.pdf"%(label, region), format='pdf', dpi=300)
plt.close()


plt.plot(np.array(times), meancorrectedMeanVolumes, 'ro', np.array(times), meancorrectedMeanVolumes)
plt.plot(np.array(times), meancorrectedMeanVolumes + meancorrectedStdVolumes, 'r--')
plt.plot(np.array(times), meancorrectedMeanVolumes - meancorrectedStdVolumes, 'r--')
plt.xlabel('time (hrs)')
plt.ylabel('cell mean/std volume after normalization by mean')
plt.title(region)
plt.savefig(savePath + "%s_cell_mean_volume_normalized_%s.pdf"%(label, region), format='pdf', dpi=300)
plt.close()

average = np.mean(meanVolumes)
plt.plot(np.array(times), meanVolumes/average-1.0)
plt.xlabel('time (hrs)')
plt.ylabel('stretch corrected: (mean volume -overall mean)/overall mean ')
plt.title(region)
plt.savefig(savePath + "%s_stretch_corrected_mean_volume_relative_error_%s.pdf"%(label, region), format='pdf', dpi=300)
plt.close()




