import cPickle
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from scipy.stats.stats import skew
import all_plant_variables as apv 

# # # # # # # # # # # # # CELLS TO INCLUDE IN STATS: CENTRAL ZONE ONLY OR CLEANED FULL LINEAGES 
radius = '30'
region = 'initial central L1' 

# split data according to: "all_data", "early_times", "late_times", "inner_zone", "outer_zone", "small_volumes", "large_volumes"
zone ="all_data"  # all_data, early_times, late_times, inner_zone, outer_zone
resultsPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_pooled_data/" + zone + "/"

# # # # # # # # # # # # # VARIABLE LABELS 
cell_vars = ["V", "V norm", "A", "Aop",  "Aip", "Aa", "cyto V"]

# # # # # # # # # # # # # FIGURES AND FIGURE VARIABLES
plantFontSize = 7
plt.locator_params(axis = 'x', nbins=4); plt.locator_params(axis = 'y', nbins=5)
plt.rcParams['xtick.labelsize'] = 20; plt.rcParams['ytick.labelsize'] = 20

# # # # # # # # # # # # # FORMAT FOR EXPORTED DATA
varBVsvarD_fits = [[" var Bir. vs var Div.", "slope", "intercept", "corr. coeff.", "p-value", " var Bir. vs Inc.", "slope", "intercept", "corr. coeff.", "p-value", "mean var birth", "coeff. var. var birth", "mean var div.", "std var division/mean var birth", "std (var S1 - var S2)/ 2(var S1 + var S2)", "med. (S_d/S_b)"]]

all_plants, all_cellCycleTimes, all_distances, all_times, all_sis_ratios_birth, all_birth_vol, all_div_vol = [], [], [], [], [], [], []
all_vars_birth, all_vars_div, all_vars_asym_div = [[] for var in cell_vars], [[] for var in cell_vars], [[] for var in cell_vars]
meanCycleTime = []

# # # # # # # # # # # # # FUNCTION THAT SPLITS DATA
def include_sample(X):
    if(zone=="inner_zone_15micron"):
        if(dist_centre[X] < 15): return True
        else: return False
    elif(zone=="outer_zone_15micron"):
        if(dist_centre[X] >= 15): return True
        else: return False
    elif(zone=="early_times"):
        if(time_at_birth[X] < np.median(time_at_birth)): return True
        else: return False
    elif(zone=="first_quarter_times"):
        if(time_at_birth[X] < np.percentile(time_at_birth, 25)): return True
        else: return False
    elif(zone=="late_times"):
        if(time_at_birth[X] >= np.median(time_at_birth)): return True
        else: return False
    elif(zone=="small_volumes"):
        if(size_vars_birth[0][X] < np.median(size_vars_birth[0])): return True
        else: return False
    elif(zone=="large_volumes"):
        if(size_vars_birth[0][X] >= np.median(size_vars_birth[0])): return True
        else: return False
    elif(zone=="large_sisters"):
        if(sis_ratio_birth[X] >= 0.0): return True
        else: return False
    elif(zone=="small_sisters"):
        if(sis_ratio_birth[X] < 0.0): return True
        else: return False
    elif(zone =="symmetric_division"):
        if(sis_ratio_birth[X]*sis_ratio_birth[X] < 0.06*0.06): return True
        else: return False
    elif(zone == "small_num_L1_neighbours_birth"):
        if(number_L1_neighbors_birth[X] < 6): return True
        else: return False
    elif(zone == "large_num_L1_neighbours_birth"):
        if(number_L1_neighbors_birth[X] > 6): return True
        else: return False
    elif(zone == "six_L1_neighbours_birth"):
        if(number_L1_neighbors_birth[X] == 6): return True
        else: return False
    elif(zone == "small_num_L1_neighbours_division"):
        if(number_L1_neighbors_division[X] < 7): return True
        else: return False
    elif(zone == "large_num_L1_neighbours_division"):
        if(number_L1_neighbors_division[X] > 7): return True
        else: return False
    elif(zone == "seven_L1_neighbours_division"):
        if(number_L1_neighbors_division[X] == 7): return True
        else: return False
    else: return True

def round_sig(x, sig=2):
    return np.round(x, sig-int(np.floor(np.log10(x)))-1)

# # # # # # # # # # # # # MAIN LOOP
for plant in range(len(apv.allPlants)):
    meanCycleTime.append([])
    rootPath = apv.allRootPaths[plant]; timesList = apv.allTimesList[plant]; ALTDataPathList = apv.allALTDataPathList[plant]; segDataPathList = apv.allSegDataPathList[plant]
    meanVolumes = apv.allMeanVolumes[plant]; L1centre = apv.allL1centre[plant]
    dataFolder = rootPath 
    times = [timesList[i][4] for i in range(0, len(timesList))]

    data_t0 = np.loadtxt(dataFolder + "Div_t0_new_" + radius + "_" + region + ".dat").astype(int) #(j0, m0, d0) , cell divides between time[j0] with mother m0 and time[j0+1] with daughter d0
    data_t1 = np.loadtxt(dataFolder + "Div_t1_new_" + radius + "_" + region + ".dat").astype(int) #(j1, m1, d1) , cell divides between time[j1] with mother m1 and time[j1+1] with daughter d1


    size_vars = [[cell_vars[i], []] for i in range(len(cell_vars))]
    size_vars_div, asymmetry_div = [[] for i in range(len(cell_vars))], [[] for i in range(len(cell_vars))]
    dist_tbirth, dist_centre = [], []
    time_since_birth, distance_over_cycle = [], []
    time_at_birth, time_of_day_at_birth  = [], []
    cell_cycle_time = []
    sis_ratio_birth = []
    number_L1_neighbors_birth = []
    number_L1_neighbors_division = []

    # for each division pair, create vectors of time, volume, normed volume, total surface area and periclinal S.A. increase between divisions 
    for i in range(0, len(data_t0)):
        print i, " from ", len(data_t0)
      
        d0 = data_t0[i][2]  # daughter cell at birth at time[j0 + 1]  
        for var in size_vars: var[1].append([])
     
        distance_over_cycle.append([]); time_since_birth.append([]); dist_evolution = []
        for j in range(data_t0[i][0] + 1, data_t1[i][0] + 1): # from time of birth to time of division
            fobj = file(segDataPathList[j])
            dataSeg = cPickle.load(fobj)
            fobj.close()
            # compute distance from L1 centre at time of birth
            if j == data_t0[i][0] + 1: 
                xSq = (dataSeg[0]['barycenter'][d0][0] - L1centre[j][0]) * (dataSeg[0]['barycenter'][d0][0] - L1centre[j][0])
                ySq = (dataSeg[0]['barycenter'][d0][1] - L1centre[j][1]) * (dataSeg[0]['barycenter'][d0][1] - L1centre[j][1])
                zSq = (dataSeg[0]['barycenter'][d0][2] - L1centre[j][2]) * (dataSeg[0]['barycenter'][d0][2] - L1centre[j][2])
                dist_tbirth.append(np.sqrt(xSq + ySq + zSq))
                time_at_birth.append(times[j])
                time_of_day_at_birth.append(times[j]%24)
                cell_cycle_time.append(times[data_t1[i][0] + 1] - times[j])
                 # compute ratio of cell volume: sister volume at birth
                fobj = file(ALTDataPathList[j-1])
                dataLineage, scores = cPickle.load(fobj)
                fobj.close()
                mother = [x[1] for x in dataLineage if x[0] == d0] 
 
                if(len(mother) > 1): print "problem with mother cell..."
                allDaughters = [x[0] for x in dataLineage if x[1] == mother[0]]
                sister = [x for x in allDaughters if x != d0]
                if(len(sister)>1): print "problem with sister cell..."

                sis_ratio_birth.append((dataSeg[0]['volumes'][d0] - dataSeg[0]['volumes'][sister[0]])/(dataSeg[0]['volumes'][d0] + dataSeg[0]['volumes'][sister[0]])) 
                number_L1_neighbors_birth.append(len([x for x in dataSeg[0]['neigbourhood'][d0] if 1 in dataSeg[0]['neigbourhood'][x]]))
            if j == data_t1[i][0]: 
                number_L1_neighbors_division.append(len([x for x in dataSeg[0]['neigbourhood'][d0] if 1 in dataSeg[0]['neigbourhood'][x]]))


            size_vars[0][1][-1].append(dataSeg[0]['volumes'][d0]); size_vars[1][1][-1].append(dataSeg[0]['volumes'][d0]/meanVolumes[j])  
            time_since_birth[-1].append(times[j] - times[data_t0[i][0] + 1])

            if(time_since_birth[-1][-1]< 0): 
                print "BUG  ... ", j, data_t0[i], time_since_birth[-1] 

            xSq = (dataSeg[0]['barycenter'][d0][0] - L1centre[j][0]) * (dataSeg[0]['barycenter'][d0][0] - L1centre[j][0])
            ySq = (dataSeg[0]['barycenter'][d0][1] - L1centre[j][1]) * (dataSeg[0]['barycenter'][d0][1] - L1centre[j][1])
            zSq = (dataSeg[0]['barycenter'][d0][2] - L1centre[j][2]) * (dataSeg[0]['barycenter'][d0][2] - L1centre[j][2])
            dist_evolution.append(np.sqrt(xSq + ySq + zSq))
            distance_over_cycle[-1].append(np.sqrt(xSq + ySq + zSq)) 

            pArea = 0; sArea = 0; oPArea = 0
            if 1 not in dataSeg[0]['neigbourhood'][d0]: print "No periclinal wall !!!!!!!!!!!!!!!!!!!!!! "
            #else: pArea = dataSeg[0]['wall_surface'][(1, d0)]
            else: pArea = dataSeg[0]['wall_surface'][(d0, 1)]

            for neighbourIndex in dataSeg[0]['neigbourhood'][d0]:
                if(neighbourIndex > d0): sArea = sArea + dataSeg[0]['wall_surface'][(neighbourIndex, d0)]
                else: sArea = sArea + dataSeg[0]['wall_surface'][(d0, neighbourIndex)]
                if 1 not in dataSeg[0]['neigbourhood'][neighbourIndex] and neighbourIndex != 1: 
                    if(neighbourIndex > d0): oPArea = oPArea + dataSeg[0]['wall_surface'][(neighbourIndex, d0)]
                    else : oPArea = oPArea + dataSeg[0]['wall_surface'][(d0, neighbourIndex)]
            size_vars[2][1][-1].append(sArea); size_vars[3][1][-1].append(pArea); size_vars[4][1][-1].append(oPArea); size_vars[5][1][-1].append(sArea - pArea - oPArea)   
            size_vars[6][1][-1].append((dataSeg[0]['volumes'][d0] - 0.06 *sArea)/meanVolumes[j])

        #find tracked d0 at next timepoint and set equal to d0. Should be only one mother because of ways cells have been selected ... 
            fobj = file(ALTDataPathList[j])
            dataLineage, scores = cPickle.load(fobj)
            fobj.close()
            newDaughters = [x[0] for x in dataLineage if x[1] == d0]
            if (j != data_t1[i][0] ) and (len(newDaughters) != 1) : print " More or less than one offspring ...  "
            if len(newDaughters) == 0 : print " Zero offspring ...  "; break
            else : d0 = newDaughters[0] 
        dist_centre.append(np.mean(dist_evolution))
        fobj = file(segDataPathList[data_t1[i][0] + 1])
        dataSeg = cPickle.load(fobj)

         
        vol, volNorm, sArea, pArea, oPArea, cytoVol = [0 for j in range(len(newDaughters))], [0 for j in range(len(newDaughters))], [0 for j in range(len(newDaughters))], [0 for j in range(len(newDaughters))], [0 for j in range(len(newDaughters))], [0 for j in range(len(newDaughters))]

        for j in range(len(newDaughters)): 
            vol[j] = dataSeg[0]['volumes'][newDaughters[j]]; volNorm[j] = dataSeg[0]['volumes'][newDaughters[j]]/meanVolumes[data_t1[i][0] + 1]
            if 1 not in dataSeg[0]['neigbourhood'][newDaughters[j]] : print "No periclinal wall !!!!!!!!!!!!!!!!!!!!!! "
            else: pArea[j] = dataSeg[0]['wall_surface'][(newDaughters[j], 1)]
            for neighbourIndex in dataSeg[0]['neigbourhood'][newDaughters[j]]:
                if(neighbourIndex > newDaughters[j]): sArea[j] = sArea[j] + dataSeg[0]['wall_surface'][(neighbourIndex, newDaughters[j])]
                else: sArea[j] = sArea[j] + dataSeg[0]['wall_surface'][(newDaughters[j], neighbourIndex)]
                if 1 not in dataSeg[0]['neigbourhood'][neighbourIndex] and neighbourIndex != 1: 
                    if(neighbourIndex > newDaughters[j]): oPArea[j] = oPArea[j] + dataSeg[0]['wall_surface'][(neighbourIndex, newDaughters[j])]
                    else : oPArea[j] = oPArea[j] + dataSeg[0]['wall_surface'][(newDaughters[j], neighbourIndex)]
            cytoVol[j] = (vol[j] - 0.06*sArea[j])/meanVolumes[data_t1[i][0] + 1]

        
        size_vars_div[0].append(np.sum(vol)); size_vars_div[1].append(np.sum(volNorm)); size_vars_div[2].append(np.sum(sArea)); size_vars_div[3].append(np.sum(pArea)); size_vars_div[4].append(np.sum(oPArea)); size_vars_div[5].append(np.sum(sArea) - np.sum(pArea) - np.sum(oPArea)); size_vars_div[6].append(np.sum(cytoVol))

        if len(newDaughters) == 2 : 
            random = 1 if np.random.rand() < 0.5 else -1
            asymmetry_div[0].append(0.5*(vol[0] - vol[1])*random/np.sum(vol)); asymmetry_div[1].append(0.5*(volNorm[0] - volNorm[1])*random/np.sum(volNorm));  asymmetry_div[2].append(0.5*(sArea[0] - sArea[1])*random/np.sum(sArea));  asymmetry_div[4].append(0.5* random*(oPArea[0] -oPArea[1])/np.sum(oPArea)); asymmetry_div[5].append(0.5* random*(sArea[0] - sArea[1] - (pArea[0] - pArea[1]) -  (oPArea[0] - oPArea[1]))/(np.sum(sArea) - np.sum(pArea) - np.sum(oPArea))) 
            if np.sum(pArea) != 0 : asymmetry_div[3].append(0.5* random*(pArea[0] - pArea[1])/np.sum(pArea))
            else : asymmetry_div[3].append(0)
            asymmetry_div[6].append(0.5*(cytoVol[0]-cytoVol[1])/np.sum(cytoVol))
        else: 
            for asymmetry in asymmetry_div: asymmetry.append(0)
            print " Aberant number of daughters =  ", len(newDaughters)

    for var in size_vars: var[1] = np.array(var[1])
    size_vars_birth = []
    for var in size_vars: size_vars_birth.append(np.array([var[1][i][0] for i in range(0, len(var[1]))]))


    print "min birth vol: ", np.amin(size_vars_birth[0]), "2.5% birth vol: ", np.percentile(size_vars_birth[0], 2.5), "5% birth vol: ", np.percentile(size_vars_birth[0], 5), "15% birth vol: ", np.percentile(size_vars_birth[0], 15), "mean birth vol: ", np.mean(size_vars_birth[0])
    print "max div vol: ", np.amax(size_vars_div[0]),  "85% birth vol: ", np.percentile(size_vars_div[0], 85), "95% birth vol: ", np.percentile(size_vars_div[0], 95), "97.5% birth vol: ", np.percentile(size_vars_div[0], 97.5), "mean division vol: ", np.mean(size_vars_div[0])

    include = []
    count, count1  = 0, 0 
    sel_plants = []; sel_cycle_time = []; sel_dist = []; sel_time_at_birth = []; sel_sis_ratio_birth = []; sel_vars_birth = [[] for entry in cell_vars]; sel_vars_div = [[] for entry in cell_vars]; sel_vars_asym_div = [[] for entry in cell_vars]
    for i in range(0, len(cell_cycle_time)):      
         if(include_sample(i)):
            sel_plants.append(plant); sel_cycle_time.append(cell_cycle_time[i])
            sel_dist.append(dist_centre[i]); sel_time_at_birth.append(time_at_birth[i])
            sel_sis_ratio_birth.append(sis_ratio_birth[i])
            include.append(True); count = count + 1
            for j in range(len(cell_vars)):
                sel_vars_birth[j].append(size_vars_birth[j][i])  
                sel_vars_div[j].append(size_vars_div[j][i])  
                sel_vars_asym_div[j].append(asymmetry_div[j][i])  
         else: include.append(False); count1 = count1 + 1
    print "Number of samples included /excluded: ", count, count1 


    mean_vars_birth = []
    for var_birth in sel_vars_birth: mean_vars_birth.append(np.mean(var_birth))
    mean_vars_div = []
    mean_vars_birth[1] = 1.0
    mean_vars_birth[6] = 1.0
    for var_div in sel_vars_div: mean_vars_div.append(np.mean(var_div))
    
    mean_period = np.mean(sel_cycle_time)
    print "plant ", plant, " mean T : ", mean_period, " sample size : ", len(sel_cycle_time) 

    for i in range(len(sel_cycle_time)):      
        all_plants.append(sel_plants[i]); 
        all_distances.append(sel_dist[i]/np.power(np.median(sel_vars_birth[0]), 0.3333333)); all_times.append(sel_time_at_birth[i])
        all_cellCycleTimes.append(sel_cycle_time[i]/mean_period) 
        all_sis_ratios_birth.append(sel_sis_ratio_birth[i])
        all_birth_vol.append(sel_vars_birth[0][i]) 
        all_div_vol.append(sel_vars_div[0][i]) 
        for j in range(len(cell_vars)):
            all_vars_birth[j].append(sel_vars_birth[j][i]/mean_vars_birth[j])  
            all_vars_div[j].append(sel_vars_div[j][i]/mean_vars_birth[j])  
            all_vars_asym_div[j].append(sel_vars_asym_div[j][i])  


print "______________________________ ALL DATA : "
 

print "min birth vol: ", np.amin(all_birth_vol), "2.5% birth vol: ", np.percentile(all_birth_vol, 2.5), "5% birth vol: ", np.percentile(all_birth_vol, 5), "15% birth vol: ", np.percentile(all_birth_vol, 15), "mean birth vol: ", np.mean(all_birth_vol)

print "max div vol: ", np.amax(all_div_vol),  "85% birth vol: ", np.percentile(all_div_vol, 85), "95% birth vol: ", np.percentile(all_div_vol, 95), "97.5% birth vol: ", np.percentile(all_div_vol, 97.5), "mean division vol: ", np.mean(all_div_vol)


Asplit = [ np.amin(all_vars_birth[1]), np.percentile(all_vars_birth[1], 25),  np.percentile(all_vars_birth[1], 50),  np.percentile(all_vars_birth[1], 75), np.amax(all_vars_birth[1])]
fig = plt.figure()
Abins = [[] for k in range(len(Asplit)-1)]; dABins = [[] for k in range(len(Asplit)-1)]
for k in range(len(all_vars_birth[1])):
    for l in range(len(Asplit)-1):
        if (all_vars_birth[1][k] >= Asplit[l] and all_vars_birth[1][k] < Asplit[l + 1]): 
            Abins[l].append(all_vars_birth[1][k]); dABins[l].append(np.array(all_cellCycleTimes[k]))   
x = [np.median(Abins[k]) for k in range(len(Abins))]; y = [np.median(dABins[k]) for k in range(len(dABins))]
yerr = [[y[k] - np.percentile(dABins[k], 25) for k in range(len(dABins))], [np.percentile(dABins[k], 75) - y[k] for k in range(len(dABins))]]     
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
#ax.set_aspect(0.7)
ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
ax.scatter(all_vars_birth[1], np.array(all_cellCycleTimes), s =100, color = 'c', marker ="o", edgecolors='None', alpha = 0.08) 
#ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
m, b = np.polyfit(all_vars_birth[1], np.array(all_cellCycleTimes), 1); r, p = pearsonr(all_vars_birth[1], np.array(all_cellCycleTimes))
vrange = [0.4 + 0.1*xx for xx in range(9)]
ax.plot(vrange, m * np.array(vrange) + b, linestyle = '-', color = 'r', linewidth = 1.3 )
leg = "$N$ = " + str(len(all_vars_birth[1])) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
ax.text(0.54, 0.79, leg, fontsize=11, color = 'k', transform=ax.transAxes)
ax.set_xlabel(  '$V_b/\mu_V$', fontsize= 25)
ax.set_ylabel( '$T/\mu_T$ ', fontsize= 25)
ax.set_ylim(0.15, 2.)
ax.set_yticks([0.5, 1.0, 1.5, 2.0])
ax.set_xticks([0.2, 0.6, 1.0, 1.4])
ax.set_xlim(0.1, 1.5)
ax.set_aspect('auto')
plt.gcf().subplots_adjust(bottom=0.15, left = 0.18)
fig.set_size_inches(5.75, 6.0)
plt.savefig(resultsPath + "%sBirthVsCellCycleTime_%s.pdf"%(cell_vars[1], region), format='pdf', dpi=200)
plt.close(fig)



Asplit = [ np.amin(all_vars_birth[6]), np.percentile(all_vars_birth[6], 25),  np.percentile(all_vars_birth[6], 50),  np.percentile(all_vars_birth[6], 75), np.amax(all_vars_birth[6])]
fig = plt.figure()
Abins = [[] for k in range(len(Asplit)-1)]; dABins = [[] for k in range(len(Asplit)-1)]
for k in range(len(all_vars_birth[6])):
    for l in range(len(Asplit)-1):
        if (all_vars_birth[6][k] >= Asplit[l] and all_vars_birth[6][k] < Asplit[l + 1]): 
            Abins[l].append(all_vars_birth[6][k]); dABins[l].append(np.array(all_cellCycleTimes[k]))   
x = [np.median(Abins[k]) for k in range(len(Abins))]; y = [np.median(dABins[k]) for k in range(len(dABins))]
yerr = [[y[k] - np.percentile(dABins[k], 25) for k in range(len(dABins))], [np.percentile(dABins[k], 75) - y[k] for k in range(len(dABins))]]     
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
#ax.set_aspect(0.7)
ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
ax.scatter(all_vars_birth[6], np.array(all_cellCycleTimes), s =100, color = 'c', marker ="o", edgecolors='None', alpha = 0.08) 
#ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
m, b = np.polyfit(all_vars_birth[6], np.array(all_cellCycleTimes), 1); r, p = pearsonr(all_vars_birth[6], np.array(all_cellCycleTimes))
vrange = [0.4 + 0.1*xx for xx in range(9)]
ax.plot(vrange, m * np.array(vrange) + b, linestyle = '-', color = 'r', linewidth = 1.3)
leg = "$N$ = " + str(len(all_vars_birth[1])) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
ax.text(0.54, 0.79, leg, fontsize=11, color = 'k', transform=ax.transAxes)
ax.set_xlabel(  'cytoplasmic $V_b/\mu_V$', fontsize= 25)
ax.set_ylabel( '$T/\mu_T$ ', fontsize= 25)
ax.set_ylim(0.15, 2.)
ax.set_yticks([0.5, 1.0, 1.5, 2.0])
ax.set_xticks([0.2, 0.6, 1.0, 1.4])
ax.set_xlim(0.1, 1.5)
ax.set_aspect('auto')
plt.gcf().subplots_adjust(bottom=0.15, left = 0.18)
fig.set_size_inches(5.75, 6.0)
plt.savefig(resultsPath + "%sBirthVsCellCycleTime_%s.pdf"%(cell_vars[6], region), format='pdf', dpi=200)
plt.close(fig)

Asplit = [ np.amin(all_distances), np.percentile(all_distances, 25),  np.percentile(all_distances, 50),  np.percentile(all_distances, 75), np.amax(all_distances)]
fig = plt.figure()
Abins = [[] for k in range(len(Asplit)-1)]; dABins = [[] for k in range(len(Asplit)-1)]
for k in range(len(all_vars_birth[1])):
    for l in range(len(Asplit)-1):
        if (all_distances[k] >= Asplit[l] and all_distances[k] < Asplit[l + 1]): 
            Abins[l].append(all_distances[k]); dABins[l].append(np.array(all_cellCycleTimes[k]))   
x = [np.median(Abins[k]) for k in range(len(Abins))]; y = [np.median(dABins[k]) for k in range(len(dABins))]
yerr = [[y[k] - np.percentile(dABins[k], 25) for k in range(len(dABins))], [np.percentile(dABins[k], 75) - y[k] for k in range(len(dABins))]]     
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
#ax.set_aspect(0.7)
ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
ax.scatter(all_distances, np.array(all_cellCycleTimes), s =100, color = 'c', marker ="o", edgecolors='None', alpha = 0.08) 
#ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
m, b = np.polyfit(all_distances, np.array(all_cellCycleTimes), 1); r, p = pearsonr(all_distances, np.array(all_cellCycleTimes))
vrange = [2 + xx for xx in range(7)]
ax.plot(vrange, m * np.array(vrange) + b, linestyle = '-', color = 'r', linewidth = 1.3 )
leg = "$N$ = " + str(len(all_vars_birth[1])) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
ax.text(0.54, 0.79, leg, fontsize=11, color = 'k', transform=ax.transAxes)
ax.set_xlabel(  'norm. distance from $O$', fontsize= 20)
ax.set_ylabel( 'norm. inter-division time', fontsize= 20)
ax.set_ylim(0.15, 2.)
ax.set_yticks([0.5, 1.0, 1.5, 2.0])
#ax.set_xticks([0.2, 0.6, 1.0, 1.4])
#ax.set_xlim(0.1, 1.5)
ax.set_aspect('auto')
plt.gcf().subplots_adjust(bottom=0.15, left = 0.18)
fig.set_size_inches(5.75, 6.0)
plt.savefig(resultsPath + "distVsCellCycleTime_%s.pdf"%(region), format='pdf', dpi=200)
plt.close(fig)
     
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.scatter(all_vars_birth[1], np.array(all_sis_ratios_birth), c= all_plants, marker ="o", edgecolors='None', alpha = 0.3) 
#ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
m, b = np.polyfit(all_vars_birth[1], np.array(all_sis_ratios_birth), 1); r, p = pearsonr(all_vars_birth[1], np.array(all_sis_ratios_birth))
leg = "$N$ = " + str(len(all_vars_birth[1])) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
ax.text(0.55, 0.77, leg, fontsize=8, color = 'k', transform=ax.transAxes)
ax.set_xlabel(  '$V_b/\mu_V$', fontsize= 15)
ax.set_ylabel( '$\\alpha_b$ ', fontsize= 15)
ax.set_xlim(0, 1.8)
ax.set_ylim(-1.0, 1.0)
ax.set_xticks([0.2, 0.6, 1.0, 1.4, 1.8])
ax.set_yticks([-1.0, 0.5, 0.0, 0.5, 1.0])
ax.set_aspect('auto')
plt.gcf().subplots_adjust(bottom=0.15, left = 0.15)
fig.set_size_inches(4.25, 4.0)
plt.savefig(resultsPath + "asymVBirthVs%sBirth_%s.pdf"%(cell_vars[1], region), format='pdf', dpi=200)
plt.close(fig)

t0 = np.min(all_times) -4 
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
#ax.set_aspect(0.7)
ax.scatter(np.array(all_times) - t0, all_cellCycleTimes, c= all_plants, marker ="o", edgecolors='None', s=100, alpha = 0.3) 
#ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
m, b = np.polyfit(np.array(all_times) - t0, all_cellCycleTimes, 1); r, p = pearsonr(np.array(all_times) - t0, all_cellCycleTimes)
leg = "$N$ = " + str(len(all_vars_birth[1])) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
ax.text(0.7, 0.75, leg, fontsize=13, color = 'k', transform=ax.transAxes)
ax.set_xlabel(  'time (h)', fontsize= 20)
ax.set_ylabel( 'norm. inter-division time ($T/\mu_T$)', fontsize= 20)
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = ( 0 - t0, 16 - t0, 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = ( 24 - t0, 24 + 16 - t0, 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = ( 2*24 - t0, 2*24 + 16 - t0, 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = ( 3*24 - t0, 3*24 + 16 - t0, 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest', alpha =0.1, extent = (16 - t0, 24 - t0, 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest', alpha =0.1, extent = (16 + 24 - t0, 24 + 24 - t0, 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest',alpha =0.1, extent = (16 + 2*24 - t0, 24 + 2 *24 - t0, 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest',alpha =0.1, extent = (16 + 3*24 -t0, 24 +  3*24 - t0, 0, 2*350))
plt.axvline(np.median(np.array(all_times) - t0), color='r', linestyle='dashed', linewidth=2)
plt.xlim(0, np.max(all_times) + 1)
plt.ylim(0,2.5)
ax.set_aspect('auto')
plt.gcf().subplots_adjust(bottom=0.14, left = 0.14)
fig.set_size_inches(10.0, 5.0)
plt.savefig(resultsPath + "tVsCellCycleTime_%s.pdf"%(region), format='pdf', dpi=200)
plt.close(fig)


for i in range(len(size_vars)):
    m, b = np.polyfit(all_vars_birth[i], all_vars_div[i], 1)
    r, p = pearsonr(all_vars_birth[i], all_vars_div[i])
    mI, bI = np.polyfit(all_vars_birth[i], np.array(all_vars_div[i]) - np.array(all_vars_birth[i]), 1)
    rI, pI = pearsonr(all_vars_birth[i], np.array(all_vars_div[i]) - np.array(all_vars_birth[i]))
    varBVsvarD_fits.append([cell_vars[i], m, b, r, p, cell_vars[i], mI, bI, rI, pI, np.mean(all_vars_birth[i]), np.std(all_vars_birth[i]), np.mean(all_vars_div[i]), np.std(all_vars_div[i]), np.std(all_vars_asym_div[i]), np.median(np.array(all_vars_div[i])/np.array(all_vars_birth[i]))])
    print cell_vars[i],",", m,",",  b,",", r, ",", p
np.savetxt(resultsPath + "varBVsvarD_fits_%s"%region + "n=" + str(len(all_vars_birth[0]))+ ".out", np.array(varBVsvarD_fits), delimiter="\t", fmt="%s")

xLimits = [[0.3,1.7 ],[0.3,1.7 ],[0.3,1.7 ],[0.3,1.7 ],[0.3,1.7 ],[0.3,1.7 ]]

for i in range(len(cell_vars)):
    print cell_vars[i]
    Asplit = [ np.amin(all_vars_birth[i]), np.percentile(all_vars_birth[i], 25),  np.percentile(all_vars_birth[i], 50),  np.percentile(all_vars_birth[i], 75), np.amax(all_vars_birth[i])]
    fig = plt.figure()
    Abins = [[] for k in range(len(Asplit)-1)]; dABins = [[] for k in range(len(Asplit)-1)]
    for k in range(len(all_vars_birth[i])):
        for l in range(len(Asplit)-1):
            if (all_vars_birth[i][k] >= Asplit[l] and all_vars_birth[i][k] < Asplit[l + 1]): 
                Abins[l].append(all_vars_birth[i][k]); dABins[l].append(all_vars_div[i][k])   
    x = [np.median(Abins[k]) for k in range(len(Abins))]; y = [np.median(dABins[k]) for k in range(len(dABins))]
    yerr = [[y[k] - np.percentile(dABins[k], 25) for k in range(len(dABins))], [np.percentile(dABins[k], 75) - y[k] for k in range(len(dABins))]]           
    ax = fig.add_subplot(1, 2, 1)
        #ax.set_aspect(0.7)
    ax.scatter(all_vars_birth[i], all_vars_div[i], s = 100, color= 'c', marker ="o", edgecolors='None', alpha = 0.08) 
   #ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
    m, b = np.polyfit(all_vars_birth[i], all_vars_div[i], 1); r, p = pearsonr(all_vars_birth[i], all_vars_div[i])
    vrange = [0.4 + 0.1*xx for xx in range(9)]
    ax.plot(vrange, m * np.array(vrange) + b, linestyle = '-', color = 'r',  linewidth = 1.3)
    ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
    leg = "$N$ = " + str(len(all_vars_birth[i])) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
    ax.text(0.07, 0.79, leg, fontsize=11, color = 'k', transform=ax.transAxes)
    if i == 1:
        ax.set_xlabel( '$V_b/\mu_V$ ', fontsize= 25)
        ax.set_ylabel( '$V_d/\mu_V$', fontsize= 25)
        ax.set_xticks([0.2, 0.6, 1.0, 1.4])
        ax.set_yticks([1.0, 1.4, 1.8, 2.2])
        ax.set_xlim(0.1, 1.5)
        ax.set_ylim(0.9, 2.25)
    if i == 6:
        ax.set_xlabel( 'cytoplasmic $V_b/\mu_V$ ', fontsize= 25)
        ax.set_ylabel( 'cytoplasmic $V_d/\mu_V$', fontsize= 25)
        ax.set_xticks([0.2, 0.6, 1.0, 1.4])
        ax.set_yticks([1.0, 1.4, 1.8, 2.2])
        ax.set_xlim(0.1, 1.5)
        ax.set_ylim(0.9, 2.25)
    else:
        ax.set_xlabel(  '${%s}_b$'%(cell_vars[i]), fontsize= 25)
        ax.set_ylabel( '${%s}_d$'%(cell_vars[i]), fontsize= 25)


    if i == 1: np.savetxt(resultsPath + "normVolb_vs_normVold_r_%s"%(radius) + ".out", [all_vars_birth[i], all_vars_div[i], all_plants], delimiter=",", fmt="%s")

    Abins = [[] for k in range(len(Asplit)-1)]; dABins = [[] for k in range(len(Asplit)-1)]
    for k in range(len(all_vars_birth[i])):
        for l in range(len(Asplit)-1):
            if (all_vars_birth[i][k] >= Asplit[l] and all_vars_birth[i][k] < Asplit[l + 1]): 
                Abins[l].append(all_vars_birth[i][k]); dABins[l].append(all_vars_div[i][k] -all_vars_birth[i][k])   
    x = [np.median(Abins[k]) for k in range(len(Abins))]; y = [np.median(dABins[k]) for k in range(len(dABins))]
    yerr = [[y[k] - np.percentile(dABins[k], 25) for k in range(len(dABins))], [np.percentile(dABins[k], 75) - y[k] for k in range(len(dABins))]]           
    ax = fig.add_subplot(1, 2, 2)
        #ax.set_aspect(0.7)
    ax.scatter(all_vars_birth[i], np.array(all_vars_div[i]) - np.array(all_vars_birth[i]), s = 100, color='c', marker ="o", edgecolors='None', alpha = 0.08) 
   #ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
    m, b = np.polyfit(all_vars_birth[i], np.array(all_vars_div[i]) - np.array(all_vars_birth[i]), 1); r, p = pearsonr(all_vars_birth[i], np.array(all_vars_div[i]) - np.array(all_vars_birth[i]))
    ax.plot(vrange, m * np.array(vrange) + b, linestyle = '-', color = 'r', linewidth = 1.3 )
    ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
    leg = "$N$ = " + str(len(all_vars_birth[i])) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
    ax.text(0.54, 0.79, leg, fontsize=11, color = 'k', transform=ax.transAxes)
    if i == 1:
        ax.set_xlabel( '$V_b/\mu_V$ ', fontsize= 25)
        ax.set_ylabel( '$\Delta/\mu_V$', fontsize= 25)
        ax.set_xticks([0.2, 0.6, 1.0, 1.4])
        ax.set_yticks([0.2, 0.6, 1.0, 1.4])
        ax.set_xlim(0.1, 1.5)
        ax.set_ylim(0.1, 1.5)
    elif i == 6: 
        ax.set_xlabel( 'cytoplasmic $V_b/\mu_V$ ', fontsize= 25)
        ax.set_ylabel( 'cytoplasmic $\Delta/\mu_V$', fontsize= 25)
        ax.set_xticks([0.2, 0.6, 1.0, 1.4])
        ax.set_yticks([0.2, 0.6, 1.0, 1.4])
        ax.set_xlim(0.1, 1.5)
        ax.set_ylim(0.1, 1.5)
    else:
        ax.set_xlabel( '${%s}_b$ '%(cell_vars[i]), fontsize= 25)
        ax.set_ylabel( '$\Delta$', fontsize= 25)


    # fig.tight_layout()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)
    fig.set_size_inches(13.0, 6.0)
    plt.gcf().subplots_adjust(bottom=0.15, left = 0.13)
    plt.savefig(resultsPath + "%sBirthVs%sDiv_%s.pdf"%(cell_vars[i], cell_vars[i],region), format='pdf', dpi=200)
    plt.close(fig)
    
