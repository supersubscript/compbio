import cPickle
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.stats import skew, skewtest
import statsmodels.api as sm
from scipy.stats.stats import pearsonr

##################################################
from config_2015_02_16_plant1 import rootPath, ALTDataPathList, segDataPathList, timesList, savePath
from indices_parameters_2015_02_16_plant1 import L1centre, diceThreshold
label = 'p18'
radius = 30
zMin = 8
region = 'L1'
################################################


times = [timesList[i][4] for i in range(0, len(timesList))]
timesLabel = [float(timesList[i][0][0:-3]) for i in range(0, len(timesList))]
timesLabel2 = [timesList[i][0][0:-3] for i in range(0, len(timesList))]
t_begin = 0 #time since beginning of experiment from when to include data

dataLineages, scores = [], []
segData = []
#record cell cycles at pairs [(j0, m0, d0), ...], [(j1,m1,d11, d12), ...] where cell is born between time[j0] with mother m0 and time[j0+1] with daughter d0, and divides between time[j1] and time[j1+1] with mother m1 and daughters d11, d12


def adjustFigAspect(fig,aspect=1):
    '''
    Adjust the subplot parameters so that the figure has the correct
    aspect ratio.
    '''
    xsize,ysize = fig.get_size_inches()
    minsize = min(xsize,ysize)
    xlim = .4*minsize/xsize
    ylim = .4*minsize/ysize
    if aspect < 1:
        xlim *= aspect
    else:
        ylim /= aspect
    fig.subplots_adjust(left=.5-xlim,
                        right=.5+xlim,
                        bottom=.5-ylim,
                        top=.5+ylim)

def innPeriWallArea(tInd, cellInd):
    ipArea = 0
    for neighbourIndex in segData[tInd][0]['neigbourhood'][cellInd]:
        if 1 not in segData[tInd][0]['neigbourhood'][neighbourIndex] and neighbourIndex != 1: 
            if(neighbourIndex < cellInd): ipArea = ipArea + segData[tInd][0]['wall_surface'][(neighbourIndex, cellInd)]
            else : ipArea = ipArea + segData[tInd][0]['wall_surface'][(cellInd, neighbourIndex)]
    return ipArea


def wallArea(tInd, cellInd):
    sArea = 0
    for neighbourIndex in segData[tInd][0]['neigbourhood'][cellInd]:
        if(neighbourIndex < cellInd): sArea = sArea + segData[tInd][0]['wall_surface'][(neighbourIndex, cellInd)]
        else: sArea = sArea + segData[tInd][0]['wall_surface'][(cellInd, neighbourIndex)]
    return sArea

def round_sig(x, sig=2):
    if x ==0 : return 0.0
    else: return np.round(x, sig-int(np.floor(np.log10(x)))-1)


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

motherRelGR, daughterRelGR = [], []


for i in range(0, len(segDataPathList)):

    #for each timepoint, identify each division of cells within radius of centre

    if i < len(segDataPathList) - 1: 
        #choose only indices that are tracked to next timepoint and are in L1
        indices_dup = [x[1] for x in dataLineages[i] if 1 in segData[i][0]['neigbourhood'][x[1]]]
        indices = [indices_dup[j] for j in range(len(indices_dup)) if indices_dup[j] not in indices_dup[j+1:]]  
    else: indices = [x[0] for x in dataLineages[i-1] if 1 in segData[i][0]['neigbourhood'][x[0]]]

    for index in indices: 
        divisionInDT = False
        xSq = (segData[i][0]['barycenter'][index][0] - L1centre[i][0]) * (segData[i][0]['barycenter'][index][0] - L1centre[i][0])
        ySq = (segData[i][0]['barycenter'][index][1] - L1centre[i][1]) * (segData[i][0]['barycenter'][index][1] - L1centre[i][1])
        zSq = (segData[i][0]['barycenter'][index][2] - L1centre[i][2]) * (segData[i][0]['barycenter'][index][2] - L1centre[i][2])
        if( np.sqrt(xSq + ySq + zSq) < radius and segData[i][0]['barycenter'][index][2] > zMin and times[i] > t_begin):                    
            if i < len(segDataPathList) - 1: 
                daughters = [x[0] for x in dataLineages[i] if x[1] == index]
                if(len(daughters) == 2): divisionInDT = True

        #for each division, compute all relative growth rates/<rel. growth rates>_neighbours of mother and daughter cell, append each pair to array
        if(divisionInDT):
            #compute mother rel g.r.
            mRGR = []
            dRGR = []
            for j in range(i-1, 0,-1):
                prevCells = [x[1] for x in dataLineages[j] if x[0] == index] #0 is daughter, 1 is mother
                #if previous cell divided in next interval, break
                if(len(prevCells)!=1): break
                if(1 not in segData[j][0]['neigbourhood'][prevCells[0]]): break
                if(len([x[0] for x in dataLineages[j] if x[1] == prevCells[0]]) > 1): break
                else:                     
                    relGR = (segData[j+1][0]['volumes'][index] - segData[j][0]['volumes'][prevCells[0]])*2/(segData[j+1][0]['volumes'][index] + segData[j][0]['volumes'][prevCells[0]])    
                    L1neighbours = [x for x in segData[j][0]['neigbourhood'][prevCells[0]] if 1 in segData[j][0]['neigbourhood'][x]]
                    relGRneighbours = []
                    for n0 in L1neighbours:
                        nextCells = [x[0] for x in dataLineages[j] if x[1] == n0]
                        if(len(nextCells)==1):
                            newRel = (segData[j+1][0]['volumes'][nextCells[0]] - segData[j][0]['volumes'][n0])*2/(segData[j+1][0]['volumes'][nextCells[0]] + segData[j][0]['volumes'][n0])
                            relGRneighbours.append(newRel)
                    relGRn = np.mean(relGRneighbours)  
                    mRGR.append(relGR/relGRn)  
                    index = prevCells[0]
          
            #compute daughter rel g.r.
            for d0 in daughters:
                for j in range(i+1, len(segDataPathList)-1):
                    nextCells = [x[0] for x in dataLineages[j] if x[1] == d0]
                    if(len(nextCells) > 1 or len(nextCells) == 0): break
                    if(1 not in segData[j+1][0]['neigbourhood'][nextCells[0]]): break
                    else: 
                        relGR = (segData[j+1][0]['volumes'][nextCells[0]] - segData[j][0]['volumes'][d0])*2/(segData[j+1][0]['volumes'][nextCells[0]] + segData[j][0]['volumes'][d0])    
                        L1neighbours = [x for x in segData[j][0]['neigbourhood'][d0] if 1 in segData[j][0]['neigbourhood'][x]]
                        d0 = nextCells[0] 
                        relGRneighbours = []
                        for n0 in L1neighbours:
                            nextCells = [x[0] for x in dataLineages[j] if x[1] == n0]
                            if(len(nextCells)==1):
                                newRel = (segData[j+1][0]['volumes'][nextCells[0]] - segData[j][0]['volumes'][n0])*2/(segData[j+1][0]['volumes'][nextCells[0]] + segData[j][0]['volumes'][n0])
                                relGRneighbours.append(newRel) 
                        relGRn = np.mean(relGRneighbours)  
                        dRGR.append(relGR/relGRn) 
                                            

            for mR in mRGR:
                for dR in dRGR:
                    if(dR < 2 and dR > 0 and mR < 2 and mR > 0):
                        motherRelGR.append(mR)
                        daughterRelGR.append(dR)
  
print len(daughterRelGR)
fig =plt.figure()   
ax = fig.add_subplot(1,1,1)
ax.scatter(motherRelGR, daughterRelGR, marker ="o", edgecolors='None',  alpha = 0.3)
#plt.(0.0, 4.0)
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(5, 3.5)
ax.set_xlim([0.05, 1.95]); ax.set_ylim([0.05, 1.95])
adjustFigAspect(fig,aspect=1.0)
plt.gcf().subplots_adjust(bottom=0.15, left = 0.15)
plt.ylabel('daughter $g_{rel}/g_{rel}^{neigh}$')
plt.xlabel('mother  $g_{rel}/g_{rel}^{neigh}}$')
plt.savefig(savePath + "%s_motherRelGRVsDaughterRelGR_r_%s_z_%s_%s_fixed_range.pdf"%(label,radius, zMin, region), format='pdf', dpi=300)
plt.close()

r, p = pearsonr(motherRelGR, daughterRelGR) 
leg = "$N$ = " + str(len(motherRelGR)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"


fig =plt.figure()   
ax = fig.add_subplot(1,1,1)
ax.hexbin(motherRelGR, daughterRelGR)
#plt.(0.0, 4.0)
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(5, 3.5)
ax.set_xlim([0.3, 1.7]); ax.set_ylim([0.3, 1.7])
adjustFigAspect(fig,aspect=1.0)
plt.gcf().subplots_adjust(bottom=0.18, left = 0.15)
ax.text(0.07, 0.75, leg, fontsize=13, color = 'w', transform=ax.transAxes)
plt.ylabel('daughter $g_{rel}/g_{rel}^{neigh}$', fontsize = 16)
plt.xlabel('mother  $g_{rel}/g_{rel}^{neigh}}$', fontsize = 16)
plt.savefig(savePath + "%s_hexbin_motherRelGRVsDaughterRelGR_r_%s_z_%s_%s_fixed_range.pdf"%(label, radius, zMin, region), format='pdf', dpi=300)
plt.close()



