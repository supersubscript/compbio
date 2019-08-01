import cPickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.stats import skew, skewtest
import statsmodels.api as sm

##################################################
from config_2015_02_16_plant15 import rootPath, ALTDataPathList, segDataPathList, nuclearVolPathList, timesList, savePath
from indices_parameters_2015_02_16_plant15 import L1centre, diceThreshold, meanVolumes

radius_in = '30'
region_in = 'initial central L1'

label = 'p15'
radius_out = 10
zMin = 8
region = 'L1'
################################################
plt.rcParams['xtick.labelsize'] = 18; plt.rcParams['ytick.labelsize'] = 18

times = [timesList[i][4] for i in range(0, len(timesList))]
timesLabel = [float(timesList[i][0][0:-3]) for i in range(0, len(timesList))]
timesLabel2 = [timesList[i][0][0:-3] for i in range(0, len(timesList))]

#cell birth and division times and indices 
data_t0 = np.loadtxt(rootPath + "Div_t0_new_" + radius_in + "_" + region_in + ".dat").astype(int) #(j0, m0, d0), cell divides between time[j0] with mother m0 and time[j0+1] with daughter d0
data_t1 = np.loadtxt(rootPath + "Div_t1_new_" + radius_in + "_" + region_in + ".dat").astype(int) #(j1, m1, d1), cell divides between time[j1] with mother m1 and time[j1+1] with daughter d1
#record cell cycles at pairs [(j0, m0, d0), ...], [(j1,m1,d11, d12), ...] where cell is born between time[j0] with mother m0 and time[j0+1] with daughter d0, and divides between time[j1] and time[j1+1] with mother m1 and daughters d11, d12

dataLineages, scores = [], []
segData = []
nucVolData = []
divisionFlag = []


def round_sig(x, sig=2):
    if x ==0 : return 0.0
    else: return np.round(x, sig-int(np.floor(np.log10(x)))-1)

def innPeriWallArea(tInd, cellInd):
    ipArea = 0
    for neighbourIndex in segData[tInd][0]['neigbourhood'][cellInd]:
        if 1 not in segData[tInd][0]['neigbourhood'][neighbourIndex] and neighbourIndex != 1: 
            if(neighbourIndex > cellInd): ipArea = ipArea + segData[tInd][0]['wall_surface'][(neighbourIndex, cellInd)]
            else : ipArea = ipArea + segData[tInd][0]['wall_surface'][(cellInd, neighbourIndex)]
    return ipArea


def wallArea(tInd, cellInd):
    sArea = 0
    for neighbourIndex in segData[tInd][0]['neigbourhood'][cellInd]:
        if(neighbourIndex > cellInd): sArea = sArea + segData[tInd][0]['wall_surface'][(neighbourIndex, cellInd)]
        else: sArea = sArea + segData[tInd][0]['wall_surface'][(cellInd, neighbourIndex)]
    return sArea


for i in range(len(ALTDataPathList)):
    fobj = file(ALTDataPathList[i])
    dataLineage, score = cPickle.load(fobj)
    fobj.close()
    dataLineages.append(dataLineage); scores.append(score)
    fobj = file (segDataPathList[i])    
    seg =  cPickle.load(fobj)
    fobj.close() 
    segData.append(seg)
    fobj = file (nuclearVolPathList[i])    
    nucVol =  cPickle.load(fobj)
    fobj.close() 
    nucVolData.append(nucVol)

fobj = file (segDataPathList[-1])    
seg =  cPickle.load(fobj)
fobj.close() 
segData.append(seg)
fobj = file (nuclearVolPathList[-1])    
nucVol =  cPickle.load(fobj)
fobj.close() 
nucVolData.append(nucVol)

all_cc_vols = []
all_cc_norm_vols = []
all_cc_opA = []
all_cc_nuc_vols = []
all_cc_norm_nuc_vols = []
all_cc_dist = []
all_cc_times = []

for i in range(0, len(data_t0)):
    cc_vols, cc_norm_vols, cc_opA, = [], [], []
    cc_dist, cc_times = [], []  
    mother = data_t0[i][1]
    for j in range(data_t0[i][0], data_t1[i][0]): 
        if(j == data_t0[i][0]) : cell = [data_t0[i][2]] 
        else: cell = [x[0] for x in dataLineages[j] if x[1] == mother] 
        if(len(cell) > 1): print "ERRRRRRORRRR!"
        else:  
            xSq = (segData[j+1][0]['barycenter'][cell[0]][0] - L1centre[j+1][0]) * (segData[j+1][0]['barycenter'][cell[0]][0] - L1centre[j+1][0])
            ySq = (segData[j+1][0]['barycenter'][cell[0]][1] - L1centre[j+1][1]) * (segData[j+1][0]['barycenter'][cell[0]][1] - L1centre[j+1][1])
            zSq = (segData[j+1][0]['barycenter'][cell[0]][2] - L1centre[j+1][2]) * (segData[j+1][0]['barycenter'][cell[0]][2] - L1centre[j+1][2])
            cc_dist.append(np.sqrt(xSq + ySq + zSq))
        mother = cell[0]         
        
    #select cells within Euclidean distance of O and sufficiently far from the base of the stack
    if( np.mean(cc_dist) < radius_out):  
        cc_nuc_vols, cc_norm_nuc_vols = [], []
        for j in range(data_t0[i][0], data_t1[i][0]): 
            if(j == data_t0[i][0]) : cell = [data_t0[i][2]] 
            else: cell = [x[0] for x in dataLineages[j] if x[1] == mother] 
            if(len(cell) > 1): print "ERRRRRRORRRR!"
            else: 
                cc_vols.append(segData[j+1][0]['volumes'][cell[0]])     
                cc_norm_vols.append(segData[j+1][0]['volumes'][cell[0]]/meanVolumes[j+1]) 
                cc_opA.append(segData[j+1][0]['wall_surface'][(cell[0],1)]) 
                cc_times.append(times[j+1]) 
                if (cell[0] in nucVolData[j+1].keys()):
                    cc_nuc_vols.append(nucVolData[j+1][cell[0]])
                    cc_norm_nuc_vols.append(nucVolData[j+1][cell[0]]/meanVolumes[j+1])
                else:
                    print "nuclear key missing ", j+1, cell[0]  
                    cc_nuc_vols.append(0)
                    cc_norm_nuc_vols.append(0)  
            if(j == data_t1[i][0] - 1):
                #append interdivision time and size at division 
                print data_t1[i]
                cc_vols.append(segData[j+2][0]['volumes'][data_t1[i][2]] + segData[j+2][0]['volumes'][data_t1[i][3]])
                cc_norm_vols.append((segData[data_t1[i][0] +1][0]['volumes'][data_t1[i][2]] + segData[j+2][0]['volumes'][data_t1[i][3]])/meanVolumes[j+2])
                cc_opA.append(segData[j+2][0]['wall_surface'][(data_t1[i][2],1)] + segData[j+2][0]['wall_surface'][(data_t1[i][3],1)] ) 
                cc_times.append(times[j+2])
                if (data_t1[i][2] in nucVolData[j+2].keys() and data_t1[i][3] in nucVolData[j+2].keys()):
                    cc_nuc_vols.append(nucVolData[j+2][data_t1[i][2]] + nucVolData[j+2][data_t1[i][3]])
                    cc_norm_nuc_vols.append((nucVolData[j+2][data_t1[i][2]] + nucVolData[j+2][data_t1[i][3]])/meanVolumes[j+2])
                else:
                    print "nuclear key missing upon division ", j+2, data_t1[i][2], data_t1[i][3]
                    cc_nuc_vols.append(0)
                    cc_norm_nuc_vols.append(0)  

            mother = cell[0]         
            
        all_cc_vols.append(cc_vols)
        all_cc_norm_vols.append(cc_norm_vols)
        all_cc_nuc_vols.append(cc_nuc_vols)
        all_cc_norm_nuc_vols.append(cc_norm_nuc_vols)
        all_cc_opA.append(cc_opA)
        all_cc_times.append(cc_times)
        all_cc_dist.append(cc_dist)

print len(all_cc_vols)
print len(all_cc_norm_vols)
print len(all_cc_nuc_vols)
print len(all_cc_norm_nuc_vols)
print len(all_cc_opA)
print len(all_cc_dist)
print len(all_cc_times)

np.savetxt(savePath + "birth_div_vols_r_%s_z_%s_%s"%(radius_out, zMin, region) + ".out", all_cc_vols, delimiter=",", fmt="%s")
np.savetxt(savePath + "birth_div_norm_vols_r_%s_z_%s_%s"%(radius_out, zMin, region) + ".out", all_cc_norm_vols, delimiter=",", fmt="%s")
np.savetxt(savePath + "birth_div_nuc_vols_r_%s_z_%s_%s"%(radius_out, zMin, region) + ".out", all_cc_nuc_vols, delimiter=",", fmt="%s")
np.savetxt(savePath + "birth_div_norm_nuc_vols_r_%s_z_%s_%s"%(radius_out, zMin, region) + ".out", all_cc_norm_nuc_vols, delimiter=",", fmt="%s")
np.savetxt(savePath + "birth_div_opAs_r_%s_z_%s_%s"%(radius_out, zMin, region) + ".out", all_cc_opA, delimiter=",", fmt="%s")
np.savetxt(savePath + "birth_div_dist_r_%s_z_%s_%s"%(radius_out, zMin, region) + ".out", all_cc_dist, delimiter=",", fmt="%s")
np.savetxt(savePath + "birth_div_times_r_%s_z_%s_%s"%(radius_out, zMin, region) + ".out", all_cc_times, delimiter=",", fmt="%s")

