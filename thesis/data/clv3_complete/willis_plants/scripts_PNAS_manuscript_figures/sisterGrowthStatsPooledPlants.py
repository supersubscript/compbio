import cPickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from scipy.stats.stats import skew
import all_plant_variables as apv 
import statsmodels.api as sm
from scipy.stats import mstats

# # # # # # # # # # # # # CELLS TO INCLUDE IN STATS: CENTRAL ZONE ONLY OR CLEANED FULL LINEAGES 
radius = "30"
# split data according to: "all_data", "early_times", "late_times", "inner_zone", "outer_zone", "small_volumes", "large_volumes"
zone = "all_data"  # all_data, early_times, late_times, inner_zone, outer_zone
resultsPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_pooled_sister_growth/" + zone + "/"

# # # # # # # # # # # # # VARIABLE LABELS 
cell_vars = ["V", "{V_n}", "A", "{A_{op}}",  "{A_{ip}}", "{A_{a}}"]
cell_vars_tex = ["$V$", "$V/\mu_V$", "$A$", "$A_{op}$",  "$A_{ip}$", "$A_a$"]

# # # # # # # # # # # # # FIGURES AND FIGURE VARIABLES
plantFontSize = 7
plt.locator_params(axis = 'x', nbins=4); plt.locator_params(axis = 'y', nbins=5)
plt.rcParams['xtick.labelsize'] = 13; plt.rcParams['ytick.labelsize'] = 13

all_plants, all_cellCycleTimes, all_distances, all_times, all_alpha_ns_neigh  = [], [[],[]], [[],[]], [[],[]], [[],[]]
all_vars_birth, all_vars_div, all_avRelGrowthRatesCCTimes, all_vars_asym_div, all_vars_evol = [[[] for var in cell_vars], [[] for var in cell_vars]], [[[] for var in cell_vars], [[] for var in cell_vars]], [[[] for var in cell_vars], [[] for var in cell_vars]], [[[] for var in cell_vars], [[] for var in cell_vars]], [[[] for var in cell_vars], [[] for var in cell_vars]]
all_alpha_b = [[[] for var in cell_vars], [[] for var in cell_vars]]

all_times_sisters = []; all_ratio_sisters = [[] for var in cell_vars]; all_sisters = [[[] for var in cell_vars],[[] for var in cell_vars]]; colourCode = []
all_asym_sisters = [[[] for var in cell_vars], [[] for var in cell_vars]];  all_av_rel_growth_sisters = [[[] for var in cell_vars], [[] for var in cell_vars]]; 

sp_time_of_day = np.array([0, 4, 8, 12, 16, 20, 24])
sp_cell_cycle = np.array([0,0.2, 0.4, 0.6, 0.8, 1.0])
sp_av_cell_cycle = np.array([0,0.2, 0.4, 0.6, 0.8, 1.5])
sp_cell_size_birth =  np.array([0, 0.85, 1.0, 1.15, 2.5])
sp_cell_asym_birth =  np.array([-1.0, -0.11, 0,  0.11, 1.0])
sp_cell_diff_neigh = np.array([-1.0, -0.11, 0.0, 0.11, 1.0])
sp_sum_cell_size_birth =  np.array([0, 1.7, 2.0, 2.30, 5.0])

all_growth_time_of_day = [[[] for t in range(len(sp_time_of_day) - 1)]  for var in cell_vars]
all_frac_growth_time_in_cell_cycle = [[[] for t in range(len(sp_cell_cycle) - 1)]  for var in cell_vars]
all_exp_growth_time_in_cell_cycle = [[[] for t in range(len(sp_cell_cycle) - 1)]  for var in cell_vars]
all_exp_growth_symDiv_small_time_in_cell_cycle = [[[] for t in range(len(sp_cell_cycle) - 1)]  for var in cell_vars]
all_exp_growth_symDiv_large_time_in_cell_cycle = [[[] for t in range(len(sp_cell_cycle) - 1)]  for var in cell_vars]
all_exp_growth_asymDiv_small_time_in_cell_cycle = [[[] for t in range(len(sp_cell_cycle) - 1)]  for var in cell_vars]
all_exp_growth_asymDiv_large_time_in_cell_cycle = [[[] for t in range(len(sp_cell_cycle) - 1)]  for var in cell_vars]
all_lin_growth_time_in_cell_cycle = [[[] for t in range(len(sp_cell_cycle) - 1)] for var in cell_vars]
all_exp_growth_cell_size_birth = [[[[] for k in range(2)] for t in range(len(sp_cell_size_birth) - 1)]  for var in cell_vars] #split according to asymmetry
all_exp_growth_cell_size_birth_sp_laplac_non_sisters = [[[[] for k in range(2)] for t in range(len(sp_cell_size_birth) - 1)]  for var in cell_vars] #split according to av difference between vol of cell and non-sisters
all_exp_growth_laplace_neighbours_asym_birth = [[[[] for k in range(2)] for t in range(len(sp_cell_diff_neigh) - 1)]  for var in cell_vars] #split according to sister asymmetry
all_exp_growth_cell_asym_birth_sp_size = [[[[] for k in range(2)] for t in range(len(sp_cell_asym_birth) - 1)]  for var in cell_vars] #split according to birth size
all_exp_growth_cell_asym_birth_sp_non_sisters = [[[[] for k in range(2)] for t in range(len(sp_cell_asym_birth) - 1)]  for var in cell_vars] #split according to non-sister assymetry
all_exp_growth_sum_cell_size_birth = [[[[] for k in range(2)] for t in range(len(sp_cell_size_birth) - 1)]  for var in cell_vars] #split according to birth size of mother cell

birth_size_outlier_number = 0
birth_size_mid_number = 0

rel_sis_size_outlier_number = 0
rel_sis_size_mid_number = 0

def round_sig(x, sig=2):
    return np.round(x, sig-int(np.floor(np.log10(x)))-1)

# # # # # # # # # # # # # MAIN LOOP
counter = -1
num_small_asym, num_large_asym = 0, 0
for plant in range(len(apv.allPlants)):
    counter = counter + 1
    rootPath = apv.allRootPaths[plant]; timesList = apv.allTimesList[plant]; ALTDataPathList = apv.allALTDataPathList[plant]; segDataPathList = apv.allSegDataPathList[plant]
    meanVolumes = apv.allMeanVolumes[plant]; L1centre = apv.allL1centre[plant]  
    region = 'initial central L1'

    dataFolder = rootPath
    times = [timesList[i][4] for i in range(0, len(timesList))]

    data_t0 = np.loadtxt(dataFolder + "Div_t0_new_" + radius + "_" + region + ".dat").astype(int) #(j0, m0, d0) , cell divides between time[j0] with mother m0 and time[j0+1] with daughter d0
    data_t1 = np.loadtxt(dataFolder + "Div_t1_new_" + radius + "_" + region + ".dat").astype(int) #(j1, m1, d1) , cell divides between time[j1] with mother m1 and time[j1+1] with daughter d1


    size_vars = [[[cell_vars[i], []] for i in range(len(cell_vars))],[[cell_vars[i], []] for i in range(len(cell_vars))]]
    size_vars_div, asymmetry_div = [[[] for i in range(len(cell_vars))],[[] for i in range(len(cell_vars))]], [[[] for i in range(len(cell_vars))], [[] for i in range(len(cell_vars))]]
    dist_tbirth, dist_centre = [[],[]], [[],[]]
    time_since_birth, distance_over_cycle = [[],[]], [[],[]]
    time_at_birth, time_of_day_at_birth  = [[],[]], [[],[]]
    cell_cycle_time = [[],[]]
    av_non_sis_neighbour_vol = [[],[]] 

    # for each division pair, create vectors of time, volume, normed volume, total surface area and periclinal S.A. increase between divisions 
    total = 0
    for i in range(0, len(data_t0)):
        print i, " from ", len(data_t0)
        sisInd = [i, -1]
          # daughter cell at birth at time[j0 + 1] 
        numSisters = 0
        for j in range(i+1, len(data_t0)): 
            if (data_t0[j][0] == data_t0[i][0] and data_t0[j][1] == data_t0[i][1] and data_t0[j][2] != data_t0[i][2]): 
                sisInd[1] = j; numSisters = numSisters + 1

        if(numSisters == 1): 
            for k in range(2):
                for var in size_vars[k]: var[1].append([])
                d0 = data_t0[sisInd[k]][2]; dSis = [data_t0[sisInd[(k+1)%2]][2]] #dSis is an array
                distance_over_cycle[k].append([]); time_since_birth[k].append([]); av_non_sis_neighbour_vol[k].append([]); dist_evolution = []                
                for j in range(data_t0[sisInd[k]][0] + 1, data_t1[sisInd[k]][0] + 1): # from time of birth to time of division
                    fobj = file(segDataPathList[j])
                    dataSeg = cPickle.load(fobj)
                    fobj.close()
                # compute distance from L1 centre at time of birth
                    if j == data_t0[sisInd[k]][0] + 1: 
                        xSq = (dataSeg[0]['barycenter'][d0][0] - L1centre[j][0]) * (dataSeg[0]['barycenter'][d0][0] - L1centre[j][0])
                        ySq = (dataSeg[0]['barycenter'][d0][1] - L1centre[j][1]) * (dataSeg[0]['barycenter'][d0][1] - L1centre[j][1])
                        zSq = (dataSeg[0]['barycenter'][d0][2] - L1centre[j][2]) * (dataSeg[0]['barycenter'][d0][2] - L1centre[j][2])
                        dist_tbirth[k].append(np.sqrt(xSq + ySq + zSq))
                        time_at_birth[k].append(times[j] - times[0])
                        time_of_day_at_birth[k].append(times[j]%24)
                        cell_cycle_time[k].append(times[data_t1[sisInd[k]][0] + 1] - times[j])

                    size_vars[k][0][1][-1].append(dataSeg[0]['volumes'][d0]); size_vars[k][1][1][-1].append(dataSeg[0]['volumes'][d0]/meanVolumes[j])  
                    time_since_birth[k][-1].append(times[j] - times[data_t0[sisInd[k]][0] + 1]) #does not include division time
                    sum_neigh, num_neigh = 0, 0
                    for neighbourIndex in dataSeg[0]['neigbourhood'][d0]: 
                        if(neighbourIndex != 1 and 1 in dataSeg[0]['neigbourhood'][neighbourIndex] and neighbourIndex not in dSis): 
                            sum_neigh = sum_neigh + dataSeg[0]['volumes'][neighbourIndex]
                            num_neigh = num_neigh + 1
                    av_non_sis_neighbour_vol[k][-1].append((dataSeg[0]['volumes'][d0] - sum_neigh/num_neigh)/(sum_neigh/num_neigh + dataSeg[0]['volumes'][d0]))


                    if(time_since_birth[k][-1][-1]< 0): 
                        print "BUG  ... ", j, data_t0[sisInd[k]], time_since_birth[-1] 

                    xSq = (dataSeg[0]['barycenter'][d0][0] - L1centre[j][0]) * (dataSeg[0]['barycenter'][d0][0] - L1centre[j][0])
                    ySq = (dataSeg[0]['barycenter'][d0][1] - L1centre[j][1]) * (dataSeg[0]['barycenter'][d0][1] - L1centre[j][1])
                    zSq = (dataSeg[0]['barycenter'][d0][2] - L1centre[j][2]) * (dataSeg[0]['barycenter'][d0][2] - L1centre[j][2])
                    dist_evolution.append(np.sqrt(xSq + ySq + zSq))
                    distance_over_cycle[k][-1].append(np.sqrt(xSq + ySq + zSq)) 

                    pArea = 0.0; sArea = 0.0; oPArea = 0.0
                    if 1 not in dataSeg[0]['neigbourhood'][d0]: print "No periclinal wall !!!!!!!!!!!!!!!!!!!!!! "
                    else: pArea = dataSeg[0]['wall_surface'][(d0, 1)]

                    for neighbourIndex in dataSeg[0]['neigbourhood'][d0]:
                        if(neighbourIndex > d0): sArea = sArea + dataSeg[0]['wall_surface'][(neighbourIndex, d0)]
                        else: sArea = sArea + dataSeg[0]['wall_surface'][(d0, neighbourIndex)]
                        if 1 not in dataSeg[0]['neigbourhood'][neighbourIndex] and neighbourIndex != 1: 
                            if(neighbourIndex > d0): oPArea = oPArea + dataSeg[0]['wall_surface'][(neighbourIndex, d0)]
                            else : oPArea = oPArea + dataSeg[0]['wall_surface'][(d0, neighbourIndex)]
                    size_vars[k][2][1][-1].append(sArea); size_vars[k][3][1][-1].append(pArea); size_vars[k][4][1][-1].append(oPArea); size_vars[k][5][1][-1].append(sArea - pArea - oPArea)   

            #find tracked d0 at next timepoint and set equal to d0. Should be only one mother because of ways cells have been selected ... 
                    fobj = file(ALTDataPathList[j])
                    dataLineage, scores = cPickle.load(fobj)
                    fobj.close()
                    newDaughters = [x[0] for x in dataLineage if x[1] == d0]; newSisDaughters = [x[0] for x in dataLineage if x[1] in dSis]
                    if (j != data_t1[sisInd[k]][0]) and (len(newDaughters) != 1) : print " More or less than one offspring ...  "
                    if len(newDaughters) == 0 : print " Zero offspring ...  "; break
                    else : d0 = newDaughters[0]; dSis =  newSisDaughters
                dist_centre[k].append(np.mean(dist_evolution))
                fobj = file(segDataPathList[data_t1[sisInd[k]][0] + 1])
                dataSeg = cPickle.load(fobj)
                
                vol, volNorm, sArea, pArea, oPArea = [0.0 for j in range(len(newDaughters))], [0.0 for j in range(len(newDaughters))], [0.0 for j in range(len(newDaughters))], [0.0 for j in range(len(newDaughters))], [0.0 for j in range(len(newDaughters))]

#TO HERE
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
        
                size_vars_div[k][0].append(np.sum(vol)); size_vars_div[k][1].append(np.sum(volNorm)); size_vars_div[k][2].append(np.sum(sArea)); size_vars_div[k][3].append(np.sum(pArea)); size_vars_div[k][4].append(np.sum(oPArea)); size_vars_div[k][5].append(np.sum(sArea) - np.sum(pArea) - np.sum(oPArea)) 

                if len(newDaughters) == 2 : 
                    random = 1 if np.random.rand() < 0.5 else -1
                    asymmetry_div[k][0].append(0.5*(vol[0] - vol[1])*random/np.sum(vol)); asymmetry_div[k][1].append(0.5*(volNorm[0] - volNorm[1])*random/np.sum(volNorm));  asymmetry_div[k][2].append(0.5*(sArea[0] - sArea[1])*random/np.sum(sArea));  
                    if np.sum(oPArea) != 0: asymmetry_div[k][4].append(0.5* random*(oPArea[0] -oPArea[1])/np.sum(oPArea)); 
                    else : asymmetry_div[k][4].append(0)
                    asymmetry_div[k][5].append(0.5* random*(sArea[0] - sArea[1] - (pArea[0] - pArea[1]) -  (oPArea[0] - oPArea[1]))/(np.sum(sArea) - np.sum(pArea) - np.sum(oPArea))) 
                    if np.sum(pArea) != 0 : asymmetry_div[k][3].append(0.5* random*(pArea[0] - pArea[1])/np.sum(pArea))
                    else : asymmetry_div[k][3].append(0)
                else: 
                    for asymmetry in asymmetry_div[k]: asymmetry.append(0)
                    print " Aberrant number of daughters =  ", len(newDaughters)   

                all_times[k].append(time_since_birth[k][-1])
                for j in range(len(cell_vars)): all_vars_evol[k][j].append(size_vars[k][j][1][-1]) 

    size_vars_birth = [[],[]]
    for var in size_vars[0]: size_vars_birth[0].append(np.array([var[1][i][0] for i in range(0, len(var[1]))]))
    for var in size_vars[1]: size_vars_birth[1].append(np.array([var[1][i][0] for i in range(0, len(var[1]))]))
 
    print "CHECK : "
    print len(size_vars_birth[0][0] + size_vars_birth[1][0]), len(size_vars_birth[0][0])
    
    mean_vars_birth = [np.mean(size_vars_birth[0][k].tolist() + size_vars_birth[1][k].tolist()) for k in range(len(cell_vars))]    
    mean_period = [np.mean(cell_cycle_time[0]), np.mean(cell_cycle_time[1])]        
    true_mean_period = np.mean(cell_cycle_time[0] + cell_cycle_time[1])

    #med_non_sis_neighbour_vol = np.median([np.mean(nsn) for nsn in av_non_sis_neighbour_vol[0]] + [np.mean(nsn) for nsn in av_non_sis_neighbour_vol[1]])

    med_abs_non_sis_diff = np.median([np.abs(x) for arrNsn in av_non_sis_neighbour_vol[0] for x in arrNsn] + [np.abs(x) for arrNsn in av_non_sis_neighbour_vol[1] for x in arrNsn])
    med_abs_non_sis_diff = 0.11 #This value was used for publication
    print "med_abs_non_sis_diff", med_abs_non_sis_diff
    
    print len(size_vars_birth[0][0])
    print len(size_vars_birth[0][0].tolist() + size_vars_birth[1][0].tolist())
    up_birth_vol = np.percentile(size_vars_birth[0][0].tolist() + size_vars_birth[1][0].tolist(), 75)
    low_birth_vol = np.percentile(size_vars_birth[0][0].tolist() + size_vars_birth[1][0].tolist(), 25)
    print "75% birth vol : ", up_birth_vol/mean_vars_birth[0]
    print "25% birth vol : ", low_birth_vol/mean_vars_birth[0]

    mean_growth_cc = [[[],[]] for var in range(len(cell_vars))]
    for i in range(len(cell_cycle_time[0])):
        for var in range(len(cell_vars)):
            mean_growth_cc[var][0].append(np.mean([(size_vars[0][var][1][i][j + 1] - size_vars[0][var][1][i][j])/(size_vars[0][var][1][i][j + 1] + size_vars[0][var][1][i][j])* 2 /(time_since_birth[0][i][j + 1] - time_since_birth[0][i][j]) for j in range(len(time_since_birth[0][i]) - 1) if (size_vars[0][var][1][i][j + 1] + size_vars[0][var][1][i][j]) != 0]))
            mean_growth_cc[var][1].append(np.mean([(size_vars[1][var][1][i][j + 1] - size_vars[1][var][1][i][j])/(size_vars[1][var][1][i][j + 1] + size_vars[1][var][1][i][j])* 2 /(time_since_birth[1][i][j + 1] - time_since_birth[1][i][j]) for j in range(len(time_since_birth[1][i]) - 1) if (size_vars[1][var][1][i][j + 1] + size_vars[1][var][1][i][j]) != 0]))
            if np.isnan(mean_growth_cc[var][0][-1]): 
                print [(size_vars[0][var][1][i][j + 1] + size_vars[0][var][1][i][j]) for j in range(len(time_since_birth[0][i]) - 1)]
                print [(time_since_birth[0][i][j + 1] - time_since_birth[0][i][j]) for j in range(len(time_since_birth[0][i]) - 1)]
                print [(size_vars[0][var][1][i][j + 1] - size_vars[0][var][1][i][j]) for j in range(len(time_since_birth[0][i]) - 1)]
            if np.isnan(mean_growth_cc[var][1][-1]): 
                print [(size_vars[1][var][1][i][j + 1] + size_vars[1][var][1][i][j]) for j in range(len(time_since_birth[1][i]) - 1)]
                print [(time_since_birth[1][i][j + 1] - time_since_birth[1][i][j]) for j in range(len(time_since_birth[1][i]) - 1)]
                print [(size_vars[1][var][1][i][j + 1] - size_vars[1][var][1][i][j]) for j in range(len(time_since_birth[1][i]) - 1)]


    print "true_mean_growth", len(mean_growth_cc[0][0]), len(mean_growth_cc[0][0] + mean_growth_cc[1][0])
    true_mean_growth = [np.mean(mean_growth_cc[var][0] + mean_growth_cc[var][1]) for var in range(len(cell_vars))]

    for i in range(len(cell_cycle_time[0])):     
        all_plants.append(plant); 
        for k in range(2):
            all_distances[k].append(dist_centre[k][i]/np.power(np.median(size_vars_birth[k][0]), 0.3333333)); all_cellCycleTimes[k].append(cell_cycle_time[k][i]/true_mean_period)
            
            for j in range(len(cell_vars)):
                if not(np.isnan(mean_growth_cc[j][k][i]/true_mean_growth[j])):
                    all_vars_birth[k][j].append(size_vars_birth[k][j][i]/mean_vars_birth[j])  
                    all_vars_div[k][j].append(size_vars_div[k][j][i]/mean_vars_birth[j])
                    all_avRelGrowthRatesCCTimes[k][j].append(mean_growth_cc[j][k][i]/true_mean_growth[j] *cell_cycle_time[k][i]/true_mean_period)
            for l in range(len(av_non_sis_neighbour_vol[k][i])):
                all_alpha_ns_neigh[k].append(av_non_sis_neighbour_vol[k][i][l])
                for j in range(len(cell_vars)):
                    all_alpha_b[k][j].append((size_vars_birth[k][j][i] - size_vars_birth[(k+1)%2][j][i])/(size_vars_birth[k][j][i] + size_vars_birth[(k+1)%2][j][i]))
                    

    birth_ratios = [size_vars_birth[0][0][k]/size_vars_birth[1][0][k] if size_vars_birth[0][0][k] > size_vars_birth[1][0][k] else size_vars_birth[1][0][k]/size_vars_birth[0][0][k] for k in range(len(size_vars_birth[0][0]))]
    med_birth_ratio = np.median(birth_ratios) 
    med_birth_ratio = 1.25 #this value was used for publication
    #print med_birth_ratio    

    for entry in birth_ratios:
        if entry > med_birth_ratio: num_large_asym = num_large_asym + 1
        if entry <= med_birth_ratio: num_small_asym = num_small_asym + 1
    #Find average non-neighbour sister size at birth


    for i in range(len(cell_cycle_time[0])):
        if(size_vars_birth[0][0][i] > size_vars_birth[1][0][i]):  s1 = 0; s2 = 1
        else: s1 = 1; s2 = 0  
        av_cycle_time = 0.5*(cell_cycle_time[0][i] + cell_cycle_time[1][i])
           
        for var in range(len(cell_vars)):
            asym_m = (size_vars_birth[s1][var][i] -size_vars_birth[s2][var][i])/(size_vars_birth[s1][var][i] + size_vars_birth[s2][var][i])
            #####update with new growth data from above? 
            growth1 = np.mean([(size_vars[s1][var][1][i][j + 1] - size_vars[s1][var][1][i][j])/(size_vars[s1][var][1][i][j + 1] + size_vars[s1][var][1][i][j])* 2 * true_mean_period/(time_since_birth[s1][i][j + 1] - time_since_birth[s1][i][j]) for j in range(len(time_since_birth[s1][i]) - 1) if (size_vars[s1][var][1][i][j + 1] + size_vars[s1][var][1][i][j]) != 0])
            growth2 = np.mean([(size_vars[s2][var][1][i][j + 1] - size_vars[s2][var][1][i][j])/(size_vars[s2][var][1][i][j + 1] + size_vars[s2][var][1][i][j])* 2 * true_mean_period/(time_since_birth[s2][i][j + 1] - time_since_birth[s2][i][j]) for j in range(len(time_since_birth[s2][i]) - 1) if (size_vars[s2][var][1][i][j + 1] + size_vars[s2][var][1][i][j]) != 0])

            if not(np.isnan(growth1) or np.isnan(growth2)):
                all_asym_sisters[s1][var].append(asym_m)
                all_asym_sisters[s2][var].append(-asym_m)
 
                all_av_rel_growth_sisters[s1][var].append(growth1)
                all_av_rel_growth_sisters[s2][var].append(growth2)
     
        if(cell_cycle_time[s1][i] > cell_cycle_time[s2][i]): 
            arr1 = np.array(time_since_birth[s2][i])/av_cycle_time
            all_times_sisters.append(arr1); colourCode.append(counter)
            for k in range(len(cell_vars)):
                arr2 =  np.array(size_vars[s1][k][1][i])[:len(size_vars[s2][k][1][i])]/np.array(size_vars[s2][k][1][i])
                arr3 = np.array(size_vars[s1][k][1][i])[:len(size_vars[s2][k][1][i])]
                arr4 = np.array(size_vars[s2][k][1][i])
                all_ratio_sisters[k].append(arr2); all_sisters[0][k].append(arr3); all_sisters[1][k].append(arr4)

        else: 
            arr1 = np.array(time_since_birth[s1][i])/av_cycle_time
            all_times_sisters.append(arr1); colourCode.append(counter)
            for k in range(len(cell_vars)):
                arr2 = np.array(size_vars[s1][k][1][i])/np.array(size_vars[s2][k][1][i])[:len(size_vars[s1][k][1][i])]
                arr3 = np.array(size_vars[s1][k][1][i]) 
                arr4 = np.array(size_vars[s2][k][1][i])[:len(size_vars[s1][k][1][i])]
                all_ratio_sisters[k].append(arr2); all_sisters[0][k].append(arr3); all_sisters[1][k].append(arr4)

        if (time_since_birth[0][i] < time_since_birth[1][i]): min_T = time_since_birth[0][i]
        else: min_T =  time_since_birth[1][i]
        for j in range(len(min_T)-1):
            for var in range(len(cell_vars)):
                for tdb in range(len(sp_sum_cell_size_birth) - 1):
                    sum_sis_cells_birth = size_vars_birth[0][0][i] + size_vars_birth[1][0][i]
                    if (sum_sis_cells_birth/mean_vars_birth[0] >= sp_sum_cell_size_birth[tdb] and sum_sis_cells_birth/mean_vars_birth[0]  < sp_sum_cell_size_birth[tdb + 1]) :
                        if(size_vars_birth[0][0][i]/size_vars_birth[1][0][i] < med_birth_ratio and size_vars_birth[0][0][i]/size_vars_birth[1][0][i] >= 1.0/med_birth_ratio): #symmetric div.
                            if(size_vars[0][var][1][i][j] + size_vars[0][var][1][i][j + 1] != 0 and size_vars[1][var][1][i][j] + size_vars[1][var][1][i][j + 1] != 0) :
                                all_exp_growth_sum_cell_size_birth[var][tdb][0].append((size_vars[0][var][1][i][j + 1] - size_vars[0][var][1][i][j] + size_vars[1][var][1][i][j + 1] - size_vars[1][var][1][i][j])/(size_vars[0][var][1][i][j + 1] + size_vars[0][var][1][i][j] + size_vars[1][var][1][i][j + 1] + size_vars[1][var][1][i][j])* 2 * true_mean_period/(time_since_birth[0][i][j + 1] - time_since_birth[0][i][j]))
                        if(size_vars_birth[0][0][i]/size_vars_birth[1][0][i] >= med_birth_ratio or size_vars_birth[0][0][i]/size_vars_birth[1][0][i] < 1.0/med_birth_ratio): #asymmetric div.
                            if(size_vars[0][var][1][i][j] + size_vars[0][var][1][i][j + 1] != 0 and size_vars[1][var][1][i][j] + size_vars[1][var][1][i][j + 1] != 0) :
                                all_exp_growth_sum_cell_size_birth[var][tdb][1].append((size_vars[0][var][1][i][j + 1] - size_vars[0][var][1][i][j] + size_vars[1][var][1][i][j + 1] - size_vars[1][var][1][i][j])/(size_vars[0][var][1][i][j + 1] + size_vars[0][var][1][i][j] + size_vars[1][var][1][i][j + 1] + size_vars[1][var][1][i][j])* 2 * true_mean_period/(time_since_birth[0][i][j + 1] - time_since_birth[0][i][j]))


  
        for sis in [0, 1]:
            oSis = 1 if sis ==0 else 0
            alpha_b = (size_vars_birth[sis][0][i] - size_vars_birth[oSis][0][i])/(size_vars_birth[sis][0][i] + size_vars_birth[oSis][0][i])
            av_laplac_non_sis_neighbour_vol = np.mean(av_non_sis_neighbour_vol[sis][i])

            if (np.abs(size_vars_birth[sis][0][i]/mean_vars_birth[0]  - 1) < 0.16) : birth_size_mid_number = birth_size_mid_number + 1
            else : birth_size_outlier_number = birth_size_outlier_number + 1

            for j in range(len(time_since_birth[sis][i]) - 1):
                for var in range(len(cell_vars)): 
                    for tdb in range(len(sp_cell_asym_birth) - 1):                        
                        if (alpha_b >= sp_cell_asym_birth[tdb] and alpha_b < sp_cell_asym_birth[tdb + 1]) :
                            #if(size_vars_birth[sis][0][i] < mean_vars_birth[0] ):
                            #if(size_vars_birth[sis][0][i]/mean_vars_birth[0]  < up_birth_vol and size_vars_birth[sis][0][i]/mean_vars_birth[0]  > low_birth_vol):
                            if(np.abs(size_vars_birth[sis][0][i]/mean_vars_birth[0]  - 1) < 0.16):
                                if (size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1] != 0 ) :
                                    #print (size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j]), (size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/
                                    all_exp_growth_cell_asym_birth_sp_size[var][tdb][0].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars[sis][var][1][i][j + 1] + size_vars[sis][var][1][i][j])* 2 * true_mean_period/(time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                                else:
                                    all_exp_growth_cell_asym_birth_sp_size[var][tdb][0].append(0.0)
                            else: 
                                if (size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1]  != 0 ):
                                    all_exp_growth_cell_asym_birth_sp_size[var][tdb][1].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars[sis][var][1][i][j + 1] + size_vars[sis][var][1][i][j])*2* true_mean_period/(time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                                else:
                                    all_exp_growth_cell_asym_birth_sp_size[var][tdb][1].append(0.0)
                            ############ add neighour dependence here  
                            if(np.abs(av_non_sis_neighbour_vol[sis][i][j]) < med_abs_non_sis_diff): 
                                if (size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1] != 0 ) :
                                    #print (size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j]), (size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/
                                    all_exp_growth_cell_asym_birth_sp_non_sisters[var][tdb][0].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars[sis][var][1][i][j + 1] + size_vars[sis][var][1][i][j])* 2 * true_mean_period/(time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                                else:
                                    all_exp_growth_cell_asym_birth_sp_non_sisters[var][tdb][0].append(0.0)
                            else: 
                                if (size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1]  != 0 ):
                                    all_exp_growth_cell_asym_birth_sp_non_sisters[var][tdb][1].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars[sis][var][1][i][j + 1] + size_vars[sis][var][1][i][j])*2* true_mean_period/(time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                                else:
                                    all_exp_growth_cell_asym_birth_sp_non_sisters[var][tdb][1].append(0.0)

                            ################  end here

                    for tdb in range(len(sp_cell_diff_neigh) - 1):                        
                        if (av_non_sis_neighbour_vol[sis][i][j] >= sp_cell_diff_neigh[tdb] and av_non_sis_neighbour_vol[sis][i][j] < sp_cell_diff_neigh[tdb + 1]) :
                            if(size_vars_birth[0][0][i]/size_vars_birth[1][0][i] < med_birth_ratio and size_vars_birth[0][0][i]/size_vars_birth[1][0][i] >= 1.0/med_birth_ratio):
                                if (size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1] != 0 ) :
                                    #print (size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j]), (size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/
                                    all_exp_growth_laplace_neighbours_asym_birth[var][tdb][0].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars[sis][var][1][i][j + 1] + size_vars[sis][var][1][i][j])* 2 * true_mean_period/(time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                                else:
                                    all_exp_growth_laplace_neighbours_asym_birth[var][tdb][0].append(0.0)
                            else: 
                                if (size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1]  != 0 ):
                                    all_exp_growth_laplace_neighbours_asym_birth[var][tdb][1].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars[sis][var][1][i][j + 1] + size_vars[sis][var][1][i][j])*2* true_mean_period/(time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                                else:
                                    all_exp_growth_laplace_neighbours_asym_birth[var][tdb][1].append(0.0)


                    for tdb in range(len(sp_cell_size_birth) - 1):
                        #if var == 1: mvb = np.mean(size_vars_birth[sis][var])
                       # else: mvb = mean_vars_birth[sis][var]
                        if (size_vars_birth[sis][0][i]/mean_vars_birth[0] >= sp_cell_size_birth[tdb] and size_vars_birth[sis][0][i]/mean_vars_birth[0]  < sp_cell_size_birth[tdb + 1]) :
                            #split according to asymmetry
                            if(size_vars_birth[0][0][i]/size_vars_birth[1][0][i] < med_birth_ratio and size_vars_birth[0][0][i]/size_vars_birth[1][0][i] >= 1.0/med_birth_ratio):
                                if (size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1])!= 0 :
                                   # print (size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j]), (size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/
                                    all_exp_growth_cell_size_birth[var][tdb][0].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars[sis][var][1][i][j + 1] + size_vars[sis][var][1][i][j])* 2 * true_mean_period/(time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                                else:
                                    all_exp_growth_cell_size_birth[var][tdb][0].append(0.0)
                            else: 
                                if (size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1])!= 0 :
                                    all_exp_growth_cell_size_birth[var][tdb][1].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars[sis][var][1][i][j + 1] + size_vars[sis][var][1][i][j])*2* true_mean_period/(time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                                else:
                                   all_exp_growth_cell_size_birth[var][tdb][1].append(0.0)
                            #split according to av. difference between sister and neighbours
                            if(av_laplac_non_sis_neighbour_vol  < 0):
                                if (size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1])!= 0 :
                                    all_exp_growth_cell_size_birth_sp_laplac_non_sisters[var][tdb][0].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars[sis][var][1][i][j + 1] + size_vars[sis][var][1][i][j])* 2 * true_mean_period/(time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                                else:
                                    all_exp_growth_cell_size_birth_sp_laplac_non_sisters[var][tdb][0].append(0.0)
                            else: 
                                if (size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1])!= 0 :
                                    all_exp_growth_cell_size_birth_sp_laplac_non_sisters[var][tdb][1].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars[sis][var][1][i][j + 1] + size_vars[sis][var][1][i][j])*2* true_mean_period/(time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                                else:
                                    all_exp_growth_cell_size_birth_sp_laplac_non_sisters[var][tdb][1].append(0.0)


                    for tdb in range(len(sp_time_of_day) - 1):
                        if ((time_at_birth[sis][i] + time_since_birth[sis][i][j])%24 >= sp_time_of_day[tdb] and (time_at_birth[sis][i] + time_since_birth[sis][i][j])%24 < sp_time_of_day[tdb + 1]) :
                            #all_growth_time_of_day[var][tdb].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars_div[sis][var][i] - size_vars_birth[sis][var][i])* cell_cycle_time[sis][i]/ (time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j])) 
                            if (size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1])!= 0 :
                                all_growth_time_of_day[var][tdb].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars[sis][var][1][i][j + 1] + size_vars[sis][var][1][i][j])* (cell_cycle_time[0][i] + cell_cycle_time[1][i])/(time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))

                    for tdb in range(len(sp_av_cell_cycle) - 1):                        
                        if (time_since_birth[sis][i][j]/av_cycle_time >= sp_av_cell_cycle[tdb] and time_since_birth[sis][i][j]/av_cycle_time < sp_av_cell_cycle[tdb + 1]) :
                            all_frac_growth_time_in_cell_cycle[var][tdb].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars_div[sis][var][i] - size_vars_birth[sis][var][i])*cell_cycle_time[sis][i]/ (time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                            """
                            if (size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1])!= 0 :
                                if(size_vars_birth[0][0][i]/size_vars_birth[1][0][i] < med_birth_ratio and size_vars_birth[0][0][i]/size_vars_birth[1][0][i] >= 1.0/med_birth_ratio):
                                    if(size_vars_birth[sis][0][i] > size_vars_birth[oSis][0][i] ): 
                                        all_exp_growth_symDiv_large_time_in_cell_cycle[var][tdb].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1])* 2 * true_mean_period /(time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                                    else:
                                        all_exp_growth_symDiv_small_time_in_cell_cycle[var][tdb].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1])* 2 * true_mean_period / (time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                                else: 
                                    if(size_vars_birth[sis][0][i] > size_vars_birth[oSis][0][i]): 
                                        all_exp_growth_asymDiv_large_time_in_cell_cycle[var][tdb].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1])*2 * true_mean_period/ (time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                                    else:
                                        all_exp_growth_asymDiv_small_time_in_cell_cycle[var][tdb].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1])*2 * true_mean_period/ (time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                            else:
                                if(size_vars_birth[0][0][i]/size_vars_birth[1][0][i] < med_birth_ratio and size_vars_birth[0][0][i]/size_vars_birth[1][0][i] >= 1.0/med_birth_ratio):
                                    if(size_vars_birth[sis][0][i] > size_vars_birth[oSis][0][i]):  all_exp_growth_symDiv_large_time_in_cell_cycle[var][tdb].append(0.0)
                                    else: all_exp_growth_symDiv_small_time_in_cell_cycle[var][tdb].append(0.0)
                                else: 
                                    if(size_vars_birth[sis][0][i] > size_vars_birth[oSis][0][i]): all_exp_growth_asymDiv_large_time_in_cell_cycle[var][tdb].append(0.0)
                                    else: all_exp_growth_asymDiv_small_time_in_cell_cycle[var][tdb].append(0.0)
                            """


                        if (time_since_birth[sis][i][j]/cell_cycle_time[sis][i] >= sp_cell_cycle[tdb] and time_since_birth[sis][i][j]/cell_cycle_time[sis][i] < sp_cell_cycle[tdb + 1]) :
                            #print sis, var, i, j, (size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j]), (size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1])/2.0,  av_cycle_time, (time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j])
                            if (size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1])!= 0 :
                                all_exp_growth_time_in_cell_cycle[var][tdb].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1])*2.0* true_mean_period/ (time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                                all_lin_growth_time_in_cell_cycle[var][tdb].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(mean_vars_birth[var]) * true_mean_period/ (time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                                
                            else:
                                all_exp_growth_time_in_cell_cycle[var][tdb].append(0.0)
                                all_lin_growth_time_in_cell_cycle[var][tdb].append(0.0)


                    for tdb in range(len(sp_cell_cycle) - 1):                        
                        if (time_since_birth[sis][i][j]/cell_cycle_time[sis][i] >= sp_cell_cycle[tdb] and time_since_birth[sis][i][j]/cell_cycle_time[sis][i] < sp_cell_cycle[tdb + 1]) :
                            if (size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1])!= 0 :
                                if(size_vars_birth[0][0][i]/size_vars_birth[1][0][i] < med_birth_ratio and size_vars_birth[0][0][i]/size_vars_birth[1][0][i] >= 1.0/med_birth_ratio):
                                    if(size_vars_birth[sis][0][i] > size_vars_birth[oSis][0][i] ): 
                                        all_exp_growth_symDiv_large_time_in_cell_cycle[var][tdb].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1])* 2 * true_mean_period /(time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                                    else:
                                        all_exp_growth_symDiv_small_time_in_cell_cycle[var][tdb].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1])* 2 * true_mean_period / (time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                                else: 
                                    if(size_vars_birth[sis][0][i] > size_vars_birth[oSis][0][i]): 
                                        all_exp_growth_asymDiv_large_time_in_cell_cycle[var][tdb].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1])* 2 * true_mean_period/ (time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                                    else:
                                        all_exp_growth_asymDiv_small_time_in_cell_cycle[var][tdb].append((size_vars[sis][var][1][i][j + 1] - size_vars[sis][var][1][i][j])/(size_vars[sis][var][1][i][j] + size_vars[sis][var][1][i][j + 1])* 2 * true_mean_period/ (time_since_birth[sis][i][j + 1] - time_since_birth[sis][i][j]))
                            else:
                                if(size_vars_birth[0][0][i]/size_vars_birth[1][0][i] < med_birth_ratio and size_vars_birth[0][0][i]/size_vars_birth[1][0][i] >= 1.0/med_birth_ratio):
                                    if(size_vars_birth[sis][0][i] > size_vars_birth[oSis][0][i]):  all_exp_growth_symDiv_large_time_in_cell_cycle[var][tdb].append(0.0)
                                    else: all_exp_growth_symDiv_small_time_in_cell_cycle[var][tdb].append(0.0)
                                else: 
                                    if(size_vars_birth[sis][0][i] > size_vars_birth[oSis][0][i]): all_exp_growth_asymDiv_large_time_in_cell_cycle[var][tdb].append(0.0)
                                    else: all_exp_growth_asymDiv_small_time_in_cell_cycle[var][tdb].append(0.0)


for h in range(len(cell_vars)):
    
   
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.boxplot(all_growth_time_of_day[h], positions = (sp_time_of_day[:-1] + sp_time_of_day[1:])*0.5, widths = 0.7)
    ax.set_xlabel(' time of day (sunrise at 0hrs)')
    ax.set_ylabel(' measure of growth rate ')
    ax.set_ylim([-2,3.5])
    plt.savefig(resultsPath + "times_day_vs_%s_growth_%s.pdf"%(cell_vars[h], radius), format='pdf', dpi=300)
    plt.close(fig)

    H, pval = mstats.kruskalwallis(all_growth_time_of_day[h][0], all_growth_time_of_day[h][1], all_growth_time_of_day[h][2], all_growth_time_of_day[h][3], all_growth_time_of_day[h][4])
    leg = "H, p = "  + str(round(H,2)) + ", " + str(round_sig(pval)) + "\n" 

    labs = []
    for t in range(len(all_growth_time_of_day[h])):
        labs.append(str(sp_time_of_day[t]) + '~' + str(sp_time_of_day[t+1]))

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    sm.graphics.violinplot(all_growth_time_of_day[h], labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs', 'violin_lw':0,
                                  'label_fontsize':'small',
                                  'label_rotation':30})
    ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = (0, 4.5, -5, 2*10))
    ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest', alpha =0.1, extent = (4.5, 6.6, -5, 2*10))
    ax.set_xlabel(' time of day (h)')
#   ax.set_ylabel(' $d%s/dt * T/\Delta_{%s}$ '%(cell_vars[h], cell_vars[h]))
    ax.set_ylabel(' $d%s/dt * T^{\, sis}/%s $'%(cell_vars[h], cell_vars[h]))
    ax.set_ylim([-2,3.5])
    ax.text(0.67, 0.85, leg, fontsize=8, color = 'k', transform=ax.transAxes)
    fig.set_size_inches(5.0, 4.2)
    plt.gcf().subplots_adjust(bottom=0.165, left = 0.13)
#   plt.savefig(resultsPath + "violin_times_day_vs_%s_growth_%s.pdf"%(cell_vars[h], radius), format='pdf', dpi=300)
    plt.savefig(resultsPath + "violin_times_day_vs_%s_exp_sis_growth_%s.pdf"%(cell_vars[h], radius), format='pdf', dpi=300)
    plt.close(fig)

    H, pval = mstats.kruskalwallis(all_frac_growth_time_in_cell_cycle[h][0], all_frac_growth_time_in_cell_cycle[h][1], all_frac_growth_time_in_cell_cycle[h][2], all_frac_growth_time_in_cell_cycle[h][3], all_frac_growth_time_in_cell_cycle[h][4])
    leg = "H, p = "  + str(round(H,2)) + ", " + str(round_sig(pval)) + "\n" 
    
    labs = []
    for t in range(len(all_frac_growth_time_in_cell_cycle[h])):
        labs.append(str(sp_cell_cycle[t]) + '~' + str(sp_cell_cycle[t+1]))

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    sm.graphics.violinplot(all_frac_growth_time_in_cell_cycle[h], labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs','violin_lw':0,
                                  'label_fontsize':'small',
                                  'label_rotation':30})
    ax.set_xlabel(' $t/T$')
    ax.set_ylabel(' $d%s/dt * T/\Delta_{%s} $'%(cell_vars[h], cell_vars[h]))
    ax.set_ylim([-1.2,4.0])
    #ax.text(0.67, 0.85, leg, fontsize=8, color = 'k', transform=ax.transAxes)
    fig.set_size_inches(5.0, 4.2)
    plt.gcf().subplots_adjust(bottom=0.165, left = 0.13)
    plt.savefig(resultsPath + "violin_cell_cycle_times_vs_%s_frac_growth_%s.pdf"%(cell_vars[h], radius), format='pdf', dpi=300)
    plt.close(fig)

    ["V", "{V_n}", "A", "{A_{op}}",  "{A_{ip}}", "{A_{a}}"]
    H, pval = mstats.kruskalwallis(all_exp_growth_time_in_cell_cycle[h][0], all_exp_growth_time_in_cell_cycle[h][1], all_exp_growth_time_in_cell_cycle[h][2], all_exp_growth_time_in_cell_cycle[h][3], all_exp_growth_time_in_cell_cycle[h][4])
    leg = "H, p = "  + str(round(H,2)) + ", " + str(round_sig(pval)) + "\n"    
    labs = []
    for t in range(len(all_exp_growth_time_in_cell_cycle[h])):
        labs.append('['+ str(sp_cell_cycle[t]) + ',' + str(sp_cell_cycle[t+1]) + ']')
        print "size array violin exp plot ", t, len(all_exp_growth_time_in_cell_cycle[h][t])
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    sm.graphics.violinplot(all_exp_growth_time_in_cell_cycle[h], labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs','violin_lw':0,
                                  'label_fontsize':'medium',
                                  'label_rotation':0})
    ax.set_xlabel(' $t/T$', fontsize = 25)
    if h == 1:
        ax.set_ylabel(' $g_{rel} \\times \, \mu_T \,  (V/\mu_V) $', fontsize = 25)
    else:
        ax.set_ylabel(' $g_{rel} \\times \, \mu_T \,  (%s) $'%cell_vars[h], fontsize = 25)
    ax.axhline(y= np.log(2), xmin=-1, xmax=7, c='c')
    ax.set_ylim([-0.8,2.5])
    #ax.text(0.67, 0.85, leg, fontsize=8, color = 'k', transform=ax.transAxes)
    fig.set_size_inches(6.5, 4.5)
    plt.gcf().subplots_adjust(bottom=0.25, left = 0.2)
    plt.savefig(resultsPath + "violin_cell_cycle_times_vs_%s_exp_growth_%s.pdf"%(cell_vars[h], radius), format='pdf', dpi=300)
    plt.close(fig)

    for t in range(len(all_lin_growth_time_in_cell_cycle[h])):
        print "size array violin lin plot ", t, len(all_lin_growth_time_in_cell_cycle[h][t])

    H, pval = mstats.kruskalwallis(all_lin_growth_time_in_cell_cycle[h][0], all_lin_growth_time_in_cell_cycle[h][1], all_lin_growth_time_in_cell_cycle[h][2], all_lin_growth_time_in_cell_cycle[h][3], all_lin_growth_time_in_cell_cycle[h][4])
    leg = "H, p = "  + str(round(H,2)) + ", " + str(round_sig(pval)) + "\n"  

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    sm.graphics.violinplot(all_lin_growth_time_in_cell_cycle[h], labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs','violin_lw':0,
                                  'label_fontsize':'medium',
                                  'label_rotation':0})
    ax.set_xlabel(' $t/T$', fontsize = 25)
    if h == 1:
        ax.set_ylabel(' $g_{abs} \\times \, \mu_T \, (V/\mu_V) $', fontsize = 25)
    else:
        ax.set_ylabel(' $g_{abs} \\times \, \mu_T  \, (%s) $'%cell_vars[h], fontsize = 25)
    ax.set_ylim([-0.8,3.8])
    ax.axhline(y= np.log(2), xmin=-1, xmax=7, c='c')
    #ax.text(0.67, 0.85, leg, fontsize=8, color = 'k', transform=ax.transAxes)
    fig.set_size_inches(6.0, 4.3)
    plt.gcf().subplots_adjust(bottom=0.25, left = 0.15)
    plt.savefig(resultsPath + "violin_cell_cycle_times_vs_%s_lin_growth_%s.pdf"%(cell_vars[h], radius), format='pdf', dpi=300)
    plt.close(fig)

    labs = []
    for t in range(len(all_exp_growth_symDiv_small_time_in_cell_cycle[h])):
        labs.append('['+ str(sp_cell_cycle[t]) + ',' + str(sp_cell_cycle[t+1]) + ']')
    
    for s1 in ('symDiv', 'asymDiv'):
        for s2 in ('small', 'large'):
            if s1 == 'symDiv' and s2 == 'small': arr = all_exp_growth_symDiv_small_time_in_cell_cycle[h]
            elif s1 == 'symDiv' and s2 == 'large': arr = all_exp_growth_symDiv_large_time_in_cell_cycle[h]
            elif s1 == 'asymDiv' and s2 == 'small': arr = all_exp_growth_asymDiv_small_time_in_cell_cycle[h]
            elif s1 == 'asymDiv' and s2 == 'large': arr = all_exp_growth_asymDiv_large_time_in_cell_cycle[h]

            H, pval = mstats.kruskalwallis(arr[0], arr[1], arr[2], arr[3], arr[4])
            print "size array violin exp plot ", s1, s2, len(arr[0]), len(arr[1]), len(arr[2]), len(arr[3]), len(arr[4])

            leg = "$H$, $p$ = "  + str(round(H,2)) + ", " + str(round_sig(pval)) + "\n"                   
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            sm.graphics.violinplot(arr, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs','violin_lw':0,
                                  'label_fontsize':'medium',
                                  'label_rotation':0})
            ax.set_xlabel(' $t/T$', fontsize = 23)
            if h == 1:
                ax.set_ylabel(' $g_{rel} \\times \mu_T \, (V/\mu_V) $', fontsize = 23)
            else:
                ax.set_ylabel(' $g_{rel} \\times \mu_T  \, (%s) $'%cell_vars[h], fontsize = 23)
            ax.set_ylim([-0.6,2.3])
            ax.set_yticks([0, 1, 2])
            ax.axhline(y= np.log(2), xmin=-1, xmax=7, c='c')
            #ax.text(0.67, 0.85, leg, fontsize=8, color = 'k', transform=ax.transAxes)
            fig.set_size_inches(6.0, 4.3)
            plt.gcf().subplots_adjust(bottom=0.3, left = 0.3)
            plt.savefig(resultsPath + "violin_cell_cycle_times_vs_%s_exp_growth_%s_%s_%s.pdf"%(cell_vars[h], s1, s2, radius), format='pdf', dpi=300)
            plt.close(fig)

    

#    eg = []
#    labs = []
#    for t in range(len(all_exp_growth_cell_size_birth[h])):
#        eg.append(all_exp_growth_cell_size_birth[h][t][0])
#        eg.append(all_exp_growth_cell_size_birth[h][t][1])
#        labs.append((sp_cell_size_birth[:-1] + sp_cell_size_birth[1:])[t]*0.5)
#        labs.append((sp_cell_size_birth[:-1] + sp_cell_size_birth[1:])[t]*0.5)


#    fig = plt.figure()
#    ax = fig.add_subplot(1, 1, 1)
#    sm.graphics.violinplot(eg, labels = labs, ax=ax, 
#                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs','violin_lw':0,
#                                  'label_fontsize':'small',
#                                  'label_rotation':30})
#    ax.imshow([[0,0], [1,1]], cmap=plt.cm.cool, interpolation = 'nearest', alpha =0.1, extent = (1.5, 2.5, -0.5, 20), aspect = 2.0)
#    ax.imshow([[0,0], [1,1]], cmap=plt.cm.cool, interpolation = 'nearest', alpha =0.1, extent = (3.5, 4.5, -.5, 20), aspect = 2.0)
#    ax.imshow([[0,0], [1,1]], cmap=plt.cm.cool, interpolation = 'nearest', alpha =0.1, extent = (5.5, 6.5, -.5, 20), aspect = 2.0)
#    ax.imshow([[0,0], [1,1]], cmap=plt.cm.cool, interpolation = 'nearest', alpha =0.1, extent = (7.5, 8.5, -.5, 20), aspect = 2.0)
#    ax.set_xlabel(' ${%s}_b$ '%cell_vars[h])
#    ax.set_ylabel(' $d%s/dt * \mu_T/%s$ '%(cell_vars[h], cell_vars[h]))
#    ax.set_ylim([-0.3, 2.0])
#    fig.set_size_inches(5.0, 4.2)
#    plt.gcf().subplots_adjust(bottom=0.145, left = 0.13)
#    plt.savefig(resultsPath + "violin_cell_size_birth_vs_%s_exp_growth_%s.pdf"%(cell_vars[h], radius), format='pdf', dpi=300)
#    plt.close(fig)

    ###############
    upylim = 2.2
    txtx = 0.6
    txty = 0.82
 
    ################## 

    eg = []
    labs = []
    for t in range(len(all_exp_growth_sum_cell_size_birth[h])):
        eg.append(all_exp_growth_sum_cell_size_birth[h][t][0])
        labs.append(str(sp_sum_cell_size_birth[t]) + '~' + str(sp_sum_cell_size_birth[t+1]))

    print "Kruskal Wallis H-test test for symmetric birth, for var:", cell_vars[h]
    H, pval = mstats.kruskalwallis(eg[0], eg[1], eg[2], eg[3])
    print H, pval
    leg = "H, p = "  + str(round(H,2)) + ", " + str(round_sig(pval)) + "\n"    

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    sm.graphics.violinplot(eg, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs','violin_lw':0,
                                  'label_fontsize':'small',
                                  'label_rotation':30})
    ax.imshow([[0,1], [0,1]], cmap=plt.cm.Greens, interpolation = 'bicubic', alpha =0.2, extent = (-1, 5, -0.3, upylim))
    ax.set_xlabel(' $({V^{\, sis1}}_b + {V^{\, sis2}}_b)/{\mu_V}_b$ ')
    ax.set_ylabel(' $d({%s}^{\, sis1} + {%s}^{\, sis2})/dt * \mu_T/({%s}^{\, sis1} + {%s}^{\, sis2})$ '%(cell_vars[h], cell_vars[h], cell_vars[h], cell_vars[h]))
    ax.set_ylim([-0.3, upylim])
    ax.axhline(y= np.log(2), xmin=-1, xmax=7, c='c')
    ax.text(txtx, txty, leg, fontsize=8, color = 'k', transform=ax.transAxes)
    fig.set_size_inches(5.0, 4.2)
    plt.gcf().subplots_adjust(bottom=0.15, left = 0.13)
    plt.savefig(resultsPath + "violin_sum_cell_size_birth_vs_%s_exp_growth_sum_sisters_%s_symmetric_birth.pdf"%(cell_vars[h], radius), format='pdf', dpi=300)
    plt.close(fig)


    eg = []
    labs = []
    for t in range(len(all_exp_growth_sum_cell_size_birth[h])):
        eg.append(all_exp_growth_sum_cell_size_birth[h][t][1])
        labs.append(str(sp_sum_cell_size_birth[t]) + '~' + str(sp_sum_cell_size_birth[t+1]))

    print "Kruskal Wallis H-test test for symmetric birth, for var:", cell_vars[h]
    H, pval = mstats.kruskalwallis(eg[0], eg[1], eg[2], eg[3])
    print H, pval
    leg = "H, p = "  + str(round(H,2)) + ", " + str(round_sig(pval)) + "\n"    

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    sm.graphics.violinplot(eg, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs','violin_lw':0,
                                  'label_fontsize':'small',
                                  'label_rotation':30})
    ax.imshow([[0,1], [0,1]], cmap=plt.cm.Greens, interpolation = 'bicubic', alpha =0.2, extent = (-1, 5, -0.3, upylim))
    ax.set_xlabel(' $({V^{\, sis1}}_b + {V^{\, sis2}}_b)/{\mu_V}_b$ ')
    ax.set_ylabel(' $d({%s}^{\, sis1} + {%s}^{\, sis2})/dt * \mu_T/({%s}^{\, sis1} + {%s}^{\, sis2})$ '%(cell_vars[h], cell_vars[h], cell_vars[h], cell_vars[h]))
    ax.set_ylim([-0.3, upylim])
    ax.axhline(y= np.log(2), xmin=-1, xmax=7, c='c')
    ax.text(txtx, txty, leg, fontsize=8, color = 'k', transform=ax.transAxes)
    fig.set_size_inches(5.0, 4.2)
    plt.gcf().subplots_adjust(bottom=0.15, left = 0.13)
    plt.savefig(resultsPath + "violin_sum_cell_size_birth_vs_%s_exp_growth_sum_sisters_%s_asymmetric_birth.pdf"%(cell_vars[h], radius), format='pdf', dpi=300)
    plt.close(fig)


    eg = []
    labs = []
    nums = []
    for t in range(len(all_exp_growth_cell_size_birth[h])):
        eg.append(all_exp_growth_cell_size_birth[h][t][0])
        nums.append(len(all_exp_growth_cell_size_birth[h][t][0]))
        labs.append(str(sp_cell_size_birth[t]) + '~' + str(sp_cell_size_birth[t+1]))

    print "Kruskal Wallis H-test test for symmetric birth, for var:", cell_vars[h]
    H, pval = mstats.kruskalwallis(eg[0], eg[1], eg[2], eg[3])
    print H, pval
    leg =  "$N$ = "+ np.array_str(np.array(nums)) + "\n" + "$H$, $p$ = "  + str(round(H,2)) + ", " + str(round_sig(pval)) + "\n"    

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    sm.graphics.violinplot(eg, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs','violin_lw':0,
                                  'label_fontsize':'small',
                                  'label_rotation':30})
    ax.imshow([[0,1], [0,1]], cmap=plt.cm.Greens, interpolation = 'bicubic', alpha =0.2, extent = (-1, 5, -0.3, upylim))
    ax.set_xlabel(' $V_b /{\mu_V}_b$ ', fontsize = 18)
    if h == 1:
        ax.set_ylabel(' $g_{rel} \\times \mu_T \, (V/\mu_V) $', fontsize = 18)
    else:
        ax.set_ylabel(' $g_{rel} \\times \mu_T  \, (%s) $'%cell_vars[h], fontsize = 18)
    ax.set_ylim([-0.3, upylim])
    ax.axhline(y= np.log(2), xmin=-1, xmax=7, c='c')
    ax.text(txtx, txty, leg, fontsize=10, color = 'k', transform=ax.transAxes)
    fig.set_size_inches(5.0, 4.2)
    plt.gcf().subplots_adjust(bottom=0.15, left = 0.13)
    plt.savefig(resultsPath + "violin_cell_size_birth_vs_%s_exp_growth_%s_symmetric_birth.pdf"%(cell_vars[h], radius), format='pdf', dpi=300)
    plt.close(fig)


    eg = []
    labs = []
    nums = []
    for t in range(len(all_exp_growth_cell_size_birth[h])):
        eg.append(all_exp_growth_cell_size_birth[h][t][1])
        nums.append(len(all_exp_growth_cell_size_birth[h][t][1]))
        labs.append(str(sp_cell_size_birth[t]) + '~' + str(sp_cell_size_birth[t+1]))
    print "Kruskal Wallis H-test test for asymmetric birth, for var:", cell_vars[h]
    H, pval = mstats.kruskalwallis(eg[0], eg[1], eg[2], eg[3])
    print H, pval
    leg =  "$N$ = "+ np.array_str(np.array(nums)) + "\n" + "$H$, $p$ = "  + str(round(H,2)) + ", " + str(round_sig(pval)) + "\n"    

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    sm.graphics.violinplot(eg, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs','violin_lw':0,
                                  'label_fontsize':'small',
                                  'label_rotation':30})
    ax.imshow([[0,1], [0,1]], cmap=plt.cm.Greens, interpolation = 'bicubic', alpha =0.2, extent = (-1, 5, -0.3, upylim))
    ax.set_xlabel(' $V_b/{\mu_V}_b$ ')
    if h == 1:
        ax.set_ylabel(' $g_{rel} \\times \mu_T \, (V/\mu_V) $', fontsize = 18)
    else:
        ax.set_ylabel(' $g_{rel} \\times \mu_T  \, (%s) $'%cell_vars[h], fontsize = 18)
    ax.set_ylim([-0.3, upylim])
    ax.axhline(y= np.log(2), xmin=-1, xmax=7, c='c')
    ax.text(txtx, txty, leg, fontsize=10, color = 'k', transform=ax.transAxes)
    fig.set_size_inches(5.0, 4.2)
    plt.gcf().subplots_adjust(bottom=0.15, left = 0.13)
    plt.savefig(resultsPath + "violin_cell_size_birth_vs_%s_exp_growth_%s_asymmetric_birth.pdf"%(cell_vars[h], radius), format='pdf', dpi=300)
    plt.close(fig)

    print "Kruskal Wallis H-test test for asymmetric birth var : ", cell_vars[h]
    H, pval = mstats.kruskalwallis(eg[0], eg[1], eg[2], eg[3])
    print H, pval



    eg = []
    labs = []
    nums = []
    for t in range(len(all_exp_growth_laplace_neighbours_asym_birth[h])):
        eg.append(all_exp_growth_laplace_neighbours_asym_birth[h][t][0])
        nums.append(len(all_exp_growth_laplace_neighbours_asym_birth[h][t][0]))
        labs.append(str(sp_cell_diff_neigh[t]) + '~' + str(sp_cell_diff_neigh[t+1]))

    print "Kruskal Wallis H-test test for symmetric birth, for var:", cell_vars[h]
    H, pval = mstats.kruskalwallis(eg[0], eg[1], eg[2], eg[3])
    print H, pval
    leg =  "$N$ = "+ np.array_str(np.array(nums)) + "\n" + "$H$, $p$ = "  + str(round(H,2)) + ", " + str(round_sig(pval)) + "\n"    

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    sm.graphics.violinplot(eg, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs','violin_lw':0,
                                  'label_fontsize':'small',
                                  'label_rotation':30})
    ax.imshow([[0,1], [0,1]], cmap=plt.cm.Greys, interpolation = 'bicubic', alpha =0.2, extent = (-1, 5, -0.3, upylim))
    ax.set_xlabel(' $(V -V^{ns-neigh}) /(V + V^{ns-neigh})$ ', fontsize = 18)
    if h == 1:
        ax.set_ylabel(' $g_{rel} \\times \mu_T \, (V/\mu_V) $', fontsize = 18)
    else:
        ax.set_ylabel(' $g_{rel} \\times \mu_T  \, (%s) $'%cell_vars[h], fontsize = 18)
    ax.set_ylim([-0.3, upylim])
    ax.axhline(y= np.log(2), xmin=-1, xmax=7, c='c')
    ax.text(txtx, txty, leg, fontsize=10, color = 'k', transform=ax.transAxes)
    fig.set_size_inches(5.0, 4.2)
    plt.gcf().subplots_adjust(bottom=0.15, left = 0.13)
    plt.savefig(resultsPath + "violin_laplace_neighbours_vs_%s_exp_growth_%s_symmetric_birth.pdf"%(cell_vars[h], radius), format='pdf', dpi=300)
    plt.close(fig)


    eg = []
    labs = []
    nums = [] 
    for t in range(len(all_exp_growth_laplace_neighbours_asym_birth[h])):
        eg.append(all_exp_growth_laplace_neighbours_asym_birth[h][t][1])
        nums.append(len(all_exp_growth_laplace_neighbours_asym_birth[h][t][1]))
        labs.append(str(sp_cell_diff_neigh[t]) + '~' + str(sp_cell_diff_neigh[t+1]))

    print "Kruskal Wallis H-test test for symmetric birth, for var:", cell_vars[h]
    H, pval = mstats.kruskalwallis(eg[0], eg[1], eg[2], eg[3])
    print H, pval
    leg = "$N$ = "+ np.array_str(np.array(nums)) + "\n" + "$H$, $p$ = "  + str(round(H,2)) + ", " + str(round_sig(pval)) + "\n"    

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    sm.graphics.violinplot(eg, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs','violin_lw':0,
                                  'label_fontsize':'small',
                                  'label_rotation':30})
    ax.imshow([[0,1], [0,1]], cmap=plt.cm.Greys, interpolation = 'bicubic', alpha =0.2, extent = (-1, 5, -0.3, upylim))
    ax.set_xlabel(' $(V -V^{ns-neigh}) /(V + V^{ns-neigh})$ ', fontsize = 18)
    if h == 1:
        ax.set_ylabel(' $g_{rel} \\times \mu_T \, (V/\mu_V) $', fontsize = 18)
    else:
        ax.set_ylabel(' $g_{rel} \\times \mu_T  \, (%s) $'%cell_vars[h], fontsize = 18)
    ax.set_ylim([-0.3, upylim])
    ax.axhline(y= np.log(2), xmin=-1, xmax=7, c='c')
    ax.text(txtx, txty, leg, fontsize=10, color = 'k', transform=ax.transAxes)
    fig.set_size_inches(5.0, 4.2)
    plt.gcf().subplots_adjust(bottom=0.15, left = 0.13)
    plt.savefig(resultsPath + "violin_laplace_neighbours_vs_%s_exp_growth_%s_asymmetric_birth.pdf"%(cell_vars[h], radius), format='pdf', dpi=300)
    plt.close(fig)


################################################
    eg = []
    labs = []
    nums = [] 
    for t in range(len(all_exp_growth_cell_size_birth_sp_laplac_non_sisters[h])):
        eg.append(all_exp_growth_cell_size_birth_sp_laplac_non_sisters[h][t][0])
        nums.append(len(all_exp_growth_cell_size_birth_sp_laplac_non_sisters[h][t][0]))
        labs.append(str(sp_cell_size_birth[t]) + '~' + str(sp_cell_size_birth[t+1]))
    print "Kruskal Wallis H-test test for symmetric birth, for var:", cell_vars[h]
    H, pval = mstats.kruskalwallis(eg[0], eg[1], eg[2], eg[3])
    print H, pval
    leg = "$N$ = "+ np.array_str(np.array(nums)) + "\n" + "$H$, $p$ = "  + str(round(2)) + ", " + str(round_sig(pval)) + "\n"    

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    sm.graphics.violinplot(eg, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs','violin_lw':0,
                                  'label_fontsize':'small',
                                  'label_rotation':30})
    ax.imshow([[0,1], [0,1]], cmap=plt.cm.Greens, interpolation = 'bicubic', alpha =0.2, extent = (-1, 5, -0.3, upylim))
    ax.set_xlabel(' $V_b/{\mu_V}_b$ ')
    ax.set_ylabel(' $d%s/dt \\times \mu_T/%s$ '%(cell_vars[h], cell_vars[h]))
    ax.set_ylim([-0.3, upylim])
    ax.axhline(y= np.log(2), xmin=-1, xmax=7, c='c')
    ax.text(txtx, txty, leg, fontsize=8, color = 'k', transform=ax.transAxes)
    fig.set_size_inches(5.0, 4.2)
    plt.gcf().subplots_adjust(bottom=0.15, left = 0.13)
    plt.savefig(resultsPath + "violin_cell_size_birth_vs_%s_exp_growth_%s_small_vol_diff_non_sisters.pdf"%(cell_vars[h], radius), format='pdf', dpi=300)
    plt.close(fig)

    eg = []
    labs = []
    nums = []
    for t in range(len(all_exp_growth_cell_size_birth_sp_laplac_non_sisters[h])):
        eg.append(all_exp_growth_cell_size_birth_sp_laplac_non_sisters[h][t][1])
        nums.append(len(all_exp_growth_cell_size_birth_sp_laplac_non_sisters[h][t][1]))
        labs.append(str(sp_cell_size_birth[t]) + '~' + str(sp_cell_size_birth[t+1]))
    print "Kruskal Wallis H-test test for symmetric birth, for var:", cell_vars[h]
    H, pval = mstats.kruskalwallis(eg[0], eg[1], eg[2], eg[3])
    print H, pval
    leg = "$N$ = "+ np.array_str(np.array(nums)) + "\n" + "$H$, $p$ = "  + str(round(2)) + ", " + str(round_sig(pval)) + "\n"    

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    sm.graphics.violinplot(eg, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs','violin_lw':0,
                                  'label_fontsize':'small',
                                  'label_rotation':30})
    ax.imshow([[0,1], [0,1]], cmap=plt.cm.Greens, interpolation = 'bicubic', alpha =0.2, extent = (-1, 5, -0.3, upylim))
    ax.set_xlabel(' $V_b/{\mu_V}_b$ ')
    ax.set_ylabel(' $d%s/dt \\times \mu_T/%s$ '%(cell_vars[h], cell_vars[h]))
    ax.set_ylim([-0.3, upylim])
    ax.axhline(y= np.log(2), xmin=-1, xmax=7, c='c')
    ax.text(txtx, txty, leg, fontsize=10, color = 'k', transform=ax.transAxes)
    fig.set_size_inches(5.0, 4.2)
    plt.gcf().subplots_adjust(bottom=0.15, left = 0.13)
    plt.savefig(resultsPath + "violin_cell_size_birth_vs_%s_exp_growth_%s_large_vol_diff_non_sisters.pdf"%(cell_vars[h], radius), format='pdf', dpi=300)
    plt.close(fig)

    eg = []
    labs = []
    nums = []
    for t in range(len(all_exp_growth_cell_asym_birth_sp_size[h])):
        eg.append(all_exp_growth_cell_asym_birth_sp_size[h][t][0])
        nums.append(len(all_exp_growth_cell_asym_birth_sp_size[h][t][0]))
        labs.append(str(sp_cell_asym_birth[t]) + '~' + str(sp_cell_asym_birth[t+1]))

    print "Kruskal Wallis H-test test for symmetric birth, for var:", cell_vars[h]
    H, pval = mstats.kruskalwallis(eg[0], eg[1], eg[2], eg[3])
    print H, pval
    leg =  "$N$ = "+ np.array_str(np.array(nums)) + "\n" + "$H$, $p$ = "  + str(round(H,2)) + ", " + str(round_sig(pval)) + "\n"    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    sm.graphics.violinplot(eg, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs','violin_lw':0, 
                                  'label_fontsize':'small',
                                  'label_rotation':30})
    ax.imshow([[0,1], [0,1]], cmap=plt.cm.Blues, interpolation = 'bicubic', alpha =0.2, extent = (-1, 5, -0.3, upylim))
    ax.set_xlabel(' $(V_b - {V^{\, sis}}_b)/({V}_b + {V^{\, sis}}_b)$', fontsize = 18)
    if h == 1:
        ax.set_ylabel(' $g_{rel} \\times \mu_T \, (V/\mu_V) $', fontsize = 18)
    else:
        ax.set_ylabel(' $g_{rel} \\times \mu_T  \, (%s) $'%cell_vars[h], fontsize = 18)
    ax.set_ylim([-0.3, upylim])
    ax.axhline(y= np.log(2), xmin=-1, xmax=7, c='c')
    ax.text(txtx, txty, leg, fontsize=10, color = 'k', transform=ax.transAxes)
    fig.set_size_inches(5.0, 4.2)
    plt.gcf().subplots_adjust(bottom=0.145, left = 0.13)
    plt.savefig(resultsPath + "violin_cell_asym_birth_vs_%s_exp_growth_%s_size_birth_outliers_removed.pdf"%(cell_vars[h], radius), format='pdf', dpi=300)
    plt.close(fig)

    eg = []
    labs = []
    nums = []
    for t in range(len(all_exp_growth_cell_asym_birth_sp_size[h])):
        eg.append(all_exp_growth_cell_asym_birth_sp_size[h][t][1])
        nums.append(len(all_exp_growth_cell_asym_birth_sp_size[h][t][1]))
        labs.append(str(sp_cell_asym_birth[t]) + '~' + str(sp_cell_asym_birth[t+1]))
    print "Kruskal Wallis H-test test for symmetric birth, for var:", cell_vars[h]
    H, pval = mstats.kruskalwallis(eg[0], eg[1], eg[2], eg[3])
    print H, pval
    leg = "$N$ = "+ np.array_str(np.array(nums)) + "\n" +  "$H$, $p$ = "  + str(round(H,2)) + ", " + str(round_sig(pval)) + "\n"    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    sm.graphics.violinplot(eg, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs','violin_lw':0,
                                  'label_fontsize':'small',
                                  'label_rotation':30})
    ax.imshow([[0,1], [0,1]], cmap=plt.cm.Blues, interpolation = 'bicubic', alpha =0.2, extent = (-1, 5, -0.3, upylim))
    ax.set_xlabel(' $(V_b - {V^{\, sis}}_b)/({V}_b + {V^{\, sis}}_b)$', fontsize = 18)
    if h == 1:
        ax.set_ylabel(' $g_{rel} \\times \mu_T \, (V/\mu_V) $', fontsize = 18)
    else:
        ax.set_ylabel(' $g_{rel} \\times \mu_T  \, (%s) $'%cell_vars[h], fontsize = 18)
    ax.set_ylim([-0.3, upylim])
    ax.axhline(y= np.log(2), xmin=-1, xmax=7, c='c')
    ax.text(txtx, txty, leg, fontsize=10, color = 'k', transform=ax.transAxes)
    fig.set_size_inches(5.0, 4.2)
    plt.gcf().subplots_adjust(bottom=0.145, left = 0.13)
    plt.savefig(resultsPath + "violin_cell_asym_birth_vs_%s_exp_growth_%s_size_birth_outliers.pdf"%(cell_vars[h], radius), format='pdf', dpi=300)
    plt.close(fig)
 
    eg = []
    labs = []
    nums = []
    for t in range(len(all_exp_growth_cell_asym_birth_sp_size[h])):
        eg.append(all_exp_growth_cell_asym_birth_sp_size[h][t][0] + all_exp_growth_cell_asym_birth_sp_size[h][t][1])
        nums.append(len(all_exp_growth_cell_asym_birth_sp_size[h][t][0] + all_exp_growth_cell_asym_birth_sp_size[h][t][1] ))
        labs.append(str(sp_cell_asym_birth[t]) + '~' + str(sp_cell_asym_birth[t+1]))
    print "Kruskal Wallis H-test test for symmetric birth, for var:", cell_vars[h]
    H, pval = mstats.kruskalwallis(eg[0], eg[1], eg[2], eg[3])
    print H, pval
    leg = "$N$ = "+ np.array_str(np.array(nums)) + "\n" +  "$H$, $p$ = "  + str(round(H,2)) + ", " + str(round_sig(pval)) + "\n"    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    sm.graphics.violinplot(eg, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs','violin_lw':0,
                                  'label_fontsize':'small',
                                  'label_rotation':30})
    ax.imshow([[0,1], [0,1]], cmap=plt.cm.Blues, interpolation = 'bicubic', alpha =0.2, extent = (-1, 5, -0.3, upylim))
    ax.set_xlabel(' $(V_b - {V^{\, sis}}_b)/({V}_b + {V^{\, sis}}_b)$', fontsize = 18)
    if h == 1:
        ax.set_ylabel(' $g_{rel} \\times \mu_T \, (V/\mu_V) $', fontsize = 18)
    else:
        ax.set_ylabel(' $g_{rel} \\times \mu_T  \, (%s) $'%cell_vars[h], fontsize = 18)
    ax.set_ylim([-0.3, upylim])
    ax.axhline(y= np.log(2), xmin=-1, xmax=7, c='c')
    ax.text(txtx, txty, leg, fontsize=10, color = 'k', transform=ax.transAxes)
    fig.set_size_inches(5.0, 4.2)
    plt.gcf().subplots_adjust(bottom=0.145, left = 0.13)
    plt.savefig(resultsPath + "violin_cell_asym_birth_vs_%s_exp_growth_%s_all_birth.pdf"%(cell_vars[h], radius), format='pdf', dpi=300)
    plt.close(fig)
    
    ###################################
    eg = []
    labs = []
    nums = []
    for t in range(len(all_exp_growth_cell_asym_birth_sp_non_sisters[h])):
        eg.append(all_exp_growth_cell_asym_birth_sp_non_sisters[h][t][0])
        nums.append(len(all_exp_growth_cell_asym_birth_sp_non_sisters[h][t][0]))
        labs.append(str(sp_cell_asym_birth[t]) + '~' + str(sp_cell_asym_birth[t+1]))

    print "Kruskal Wallis H-test test for symmetric birth, for var:", cell_vars[h]
    H, pval = mstats.kruskalwallis(eg[0], eg[1], eg[2], eg[3])
    print H, pval
    leg =  "$N$ = "+ np.array_str(np.array(nums)) + "\n" + "$H$, $p$ = "  + str(round(H,2)) + ", " + str(round_sig(pval)) + "\n"    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    sm.graphics.violinplot(eg, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs','violin_lw':0, 
                                  'label_fontsize':'small',
                                  'label_rotation':30})
    ax.imshow([[0,1], [0,1]], cmap=plt.cm.Blues, interpolation = 'bicubic', alpha =0.2, extent = (-1, 5, -0.3, upylim))
    ax.set_xlabel(' $(V_b - {V^{\, sis}}_b)/({V}_b + {V^{\, sis}}_b)$', fontsize = 18)
    if h == 1:
        ax.set_ylabel(' $g_{rel} \\times \mu_T \, (V/\mu_V) $', fontsize = 18)
    else:
        ax.set_ylabel(' $g_{rel} \\times \mu_T  \, (%s) $'%cell_vars[h], fontsize = 18)
    ax.set_ylim([-0.3, upylim])
    ax.axhline(y= np.log(2), xmin=-1, xmax=7, c='c')
    ax.text(txtx, txty, leg, fontsize=10, color = 'k', transform=ax.transAxes)
    fig.set_size_inches(5.0, 4.2)
    plt.gcf().subplots_adjust(bottom=0.145, left = 0.13)
    plt.savefig(resultsPath + "violin_cell_asym_birth_vs_%s_exp_growth_%s_small_asym_non_sisters.pdf"%(cell_vars[h], radius), format='pdf', dpi=300)
    plt.close(fig)


    eg = []
    labs = []
    nums = []
    for t in range(len(all_exp_growth_cell_asym_birth_sp_non_sisters[h])):
        eg.append(all_exp_growth_cell_asym_birth_sp_non_sisters[h][t][1])
        nums.append(len(all_exp_growth_cell_asym_birth_sp_non_sisters[h][t][1]))
        labs.append(str(sp_cell_asym_birth[t]) + '~' + str(sp_cell_asym_birth[t+1]))

    print "Kruskal Wallis H-test test for symmetric birth, for var:", cell_vars[h]
    H, pval = mstats.kruskalwallis(eg[0], eg[1], eg[2], eg[3])
    print H, pval
    leg =  "$N$ = "+ np.array_str(np.array(nums)) + "\n" + "$H$, $p$ = "  + str(round(H,2)) + ", " + str(round_sig(pval)) + "\n"    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    sm.graphics.violinplot(eg, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs','violin_lw':0, 
                                  'label_fontsize':'small',
                                  'label_rotation':30})
    ax.imshow([[0,1], [0,1]], cmap=plt.cm.Blues, interpolation = 'bicubic', alpha =0.2, extent = (-1, 5, -0.3, upylim))
    ax.set_xlabel(' $(V_b - {V^{\, sis}}_b)/({V}_b + {V^{\, sis}}_b)$', fontsize =18)
    if h == 1:
        ax.set_ylabel(' $g_{rel} \\times \mu_T \, (V/\mu_V) $', fontsize = 18)
    else:
        ax.set_ylabel(' $g_{rel} \\times \mu_T  \, (%s) $'%cell_vars[h], fontsize = 18)
    ax.set_ylim([-0.3, upylim])
    ax.axhline(y= np.log(2), xmin=-1, xmax=7, c='c')
    ax.text(txtx, txty, leg, fontsize=10, color = 'k', transform=ax.transAxes)
    fig.set_size_inches(5.0, 4.2)
    plt.gcf().subplots_adjust(bottom=0.145, left = 0.13)
    plt.savefig(resultsPath + "violin_cell_asym_birth_vs_%s_exp_growth_%s_large_asym_non_sisters.pdf"%(cell_vars[h], radius), format='pdf', dpi=300)
    plt.close(fig)


    ##################################################
  
    derivs = []
    deriv_ratios = [] 
    counter = 0
    for i in range(len(all_times_sisters)):
        if(len(all_times_sisters[i]) > 1):
            if(all_sisters[0][0][i][0]/all_sisters[1][0][i][0] > 1.4 ):
                counter = counter + 1
                dS1 = np.array(all_sisters[0][h][i])[1:] - np.array(all_sisters[0][h][i])[:-1]
                dS2 = np.array(all_sisters[1][h][i])[1:] - np.array(all_sisters[1][h][i])[:-1]
                derivs.append(dS1/dS2) 
                deriv_ratios.append((np.array(all_ratio_sisters[h][i])[1:] + np.array(all_ratio_sisters[h][i])[:-1])*0.5)
    smooth_derivs = [(derivs[i][:-1] + derivs[i][1:])*0.5 for i in range(len(derivs))]
    smooth_deriv_ratios = [(deriv_ratios[i][:-1] + deriv_ratios[i][1:])*0.5 for i in range(len(deriv_ratios))]
    print "Fraction of cell cycles in deriv plots:", float(counter)/float(len(all_times_sisters))

    all_deriv_ratios = [j for i in deriv_ratios for j in i]
    all_derivs = [j for i in derivs for j in i]

    Asplit = [ np.amin(all_deriv_ratios), np.percentile(all_deriv_ratios, 20), np.percentile(all_deriv_ratios, 40),  np.percentile(all_deriv_ratios, 60),  np.percentile(all_deriv_ratios, 80), np.amax(all_deriv_ratios)]
    Abins = [[] for k in range(len(Asplit)-1)]; dABins = [[] for k in range(len(Asplit)-1)]
    for k in range(len(all_deriv_ratios)):
        for l in range(len(Asplit)-1):
            if (all_deriv_ratios[k] >= Asplit[l] and all_deriv_ratios[k] < Asplit[l + 1]): 
                Abins[l].append(all_deriv_ratios[k]); dABins[l].append(all_derivs[k])   
    x = [np.log(np.median(Abins[k])) for k in range(len(Abins))]; y = [np.log(np.median(dABins[k])) for k in range(len(dABins))]
    yerr = [[y[k] - np.log(np.percentile(dABins[k], 25)) for k in range(len(dABins))], [np.log(np.percentile(dABins[k], 75)) - y[k] for k in range(len(dABins))]]           

    m, b = np.polyfit(x, y, 1); r, p = pearsonr(np.log(np.array(all_deriv_ratios)), np.log(np.array(all_derivs))) 
    #m, b = np.polyfit(np.log(np.array(all_deriv_ratios)), np.log(np.array(all_derivs)), 1); r, p = pearsonr(np.log(np.array(all_deriv_ratios)), np.log(np.array(all_derivs)))
    leg = "n = " + str(len(all_deriv_ratios)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "c.c., p = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"    

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    color = iter(plt.cm.rainbow(np.linspace(0,1,len(deriv_ratios))))
    for i in range(len(deriv_ratios)):
        c= next(color)
        ax.plot(np.log(np.array(deriv_ratios[i])), np.log(np.array(derivs[i])), c= c, alpha = 0.3)
    ax.errorbar(x, y, yerr =  yerr, linestyle='-', color = 'k', marker='o')
    #ax.plot(x, np.power(np.array(x),0.5), linestyle='-', color = 'g')
    ax.text(0.07, 0.85, leg, fontsize=6, color = 'k', transform=ax.transAxes)
    ax.set_xlabel( ' log(%s1/ %s2) '%(cell_vars[h], cell_vars[h]), fontsize= 8.0)
    ax.set_ylabel( ' log(d %s1/ d %s2) '%(cell_vars[h], cell_vars[h]), fontsize= 8.0)#ax.set_xlim([-0.2, 1.4])
    ax.set_xlim([0.0, 1.4])
    ax.set_ylim([-2.0, 4.0])
    plt.savefig(resultsPath + "%s_deriv_sisters_%s.pdf"%(cell_vars[h],radius), format='pdf', dpi=200)
    plt.close(fig)  

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    color = iter(plt.cm.rainbow(np.linspace(0,1,len(smooth_deriv_ratios))))
    for i in range(len(smooth_deriv_ratios)):
        c= next(color)
        ax.plot(smooth_deriv_ratios[i], smooth_derivs[i], c= c, alpha = 0.3)
    #ax.plot(x, 1.0*np.array(x), linestyle='-', color = 'g')
    #ax.plot(x, 0.0*np.array(x), linestyle='-', color = 'g')
    ax.set_xlabel( ' %s1/ %s2  '%(cell_vars[h], cell_vars[h]), fontsize= 8.0)
    ax.set_ylabel( 'd %s1/ d %s2  '%(cell_vars[h], cell_vars[h]), fontsize= 8.0)
    ax.set_xlim([0.5, 4.0])
    ax.set_ylim([-1.0, 5.0])
    plt.savefig(resultsPath + "%s_smooth_deriv_sisters_notlog_%s.pdf"%(cell_vars[h], radius), format='pdf', dpi=200)
    plt.close(fig)  

    ylimits = [[0.7,4.3], [0.7,4.3], [0.7, np.power(4.3, 0.666)], [0.7, np.power(4.3, 0.87)], [0.7, np.power(4.3, 1)],  [0.7, np.power(4.3, 0.5)]]
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    color = iter(plt.cm.rainbow(np.linspace(0,1,len(all_times_sisters))))
    for i in range(len(all_times_sisters)):
        c= next(color)
        ax.plot(all_times_sisters[i], all_ratio_sisters[h][i], c= 'c', alpha = 0.3) 
    ax.set_ylabel( ' $%s^{\, sisL}/%s^{\, sisS}$'%(cell_vars[h], cell_vars[h]), fontsize= 25.0)
    ax.set_xlabel( '$t/T^{\, sis}$', fontsize= 25.0)
    ax.set_ylim(ylimits[h])
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
    ax.set_yticks([1,2,3,4])
    fig.set_size_inches(7.5, 6.0)
    plt.gcf().subplots_adjust(bottom=0.2, left = 0.13)
    plt.savefig(resultsPath + "%s_sisters_%s.pdf"%(cell_vars[h],radius), format='pdf', dpi=200)
    plt.close(fig) 

    all_ratios_birth = [all_ratio_sisters[h][i][0] for i in range(len(all_ratio_sisters[h])) if all_times_sisters[i][-1] != 0 and all_ratio_sisters[h][i][-1]/all_ratio_sisters[h][i][0] != 0 and not np.isinf(np.log(all_ratio_sisters[h][i][-1]/all_ratio_sisters[h][i][0]))]
    all_slopes = [(all_ratio_sisters[h][i][-1] - all_ratio_sisters[h][i][0])/all_times_sisters[i][-1]  for i in range(len(all_ratio_sisters[h])) if all_times_sisters[i][-1] != 0 and all_ratio_sisters[h][i][-1]/all_ratio_sisters[h][i][0] != 0 and not np.isinf(np.log(all_ratio_sisters[h][i][-1]/all_ratio_sisters[h][i][0]))] 
    all_log_ratios_slopes = [np.log(all_ratio_sisters[h][i][-1]/all_ratio_sisters[h][i][0])/all_times_sisters[i][-1]  for i in range(len(all_ratio_sisters[h])) if all_times_sisters[i][-1] != 0 and all_ratio_sisters[h][i][-1]/all_ratio_sisters[h][i][0] != 0 and not np.isinf(np.log(all_ratio_sisters[h][i][-1]/all_ratio_sisters[h][i][0]))]
    all_asym_birth = [(all_ratio_sisters[h][i][0] -1)/(all_ratio_sisters[h][i][0] + 1) for i in range(len(all_ratio_sisters[h])) if all_times_sisters[i][-1] != 0 and all_ratio_sisters[h][i][-1]/all_ratio_sisters[h][i][0] != 0 and not np.isinf(np.log(all_ratio_sisters[h][i][-1]/all_ratio_sisters[h][i][0]))]
    newcolourCode = [colourCode[i] for i in range(len(all_ratio_sisters[h])) if all_times_sisters[i][-1] != 0 and all_ratio_sisters[h][i][-1]/all_ratio_sisters[h][i][0] != 0 and not np.isinf(np.log(all_ratio_sisters[h][i][-1]/all_ratio_sisters[h][i][0]))]


    Asplit = [0, 0.05, 0.11, 0.17, 1]

    Abins = [[] for k in range(len(Asplit)-1)]; dABins = [[] for k in range(len(Asplit)-1)]
    for k in range(len(all_asym_birth)):
        for l in range(len(Asplit)-1):
            if (all_asym_birth[k] >= Asplit[l] and all_asym_birth[k] < Asplit[l + 1]): 
                Abins[l].append(all_asym_birth[k]); dABins[l].append(all_log_ratios_slopes[k])   
    x = [np.median(Abins[k]) for k in range(len(Abins))]; y = [np.median(dABins[k]) for k in range(len(dABins))]
    yerr = [[y[k] - np.percentile(dABins[k], 25) for k in range(len(dABins))], [np.percentile(dABins[k], 75) - y[k] for k in range(len(dABins))]]           
    m, b = np.polyfit(x, y, 1); r, p = pearsonr(np.array(all_asym_birth), np.array(all_log_ratios_slopes)) 
    leg = "n = " + str(len(all_asym_birth)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b, 2)) + "\n"   + "c.c., p = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"    

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(all_asym_birth, all_log_ratios_slopes, color= "c",  edgecolors='None', s=65, alpha = 0.1) 
    ax.plot(x, m * np.array(x) + b, linestyle = '-', color = 'r', linewidth = 1.0 )
    ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
    ax.set_ylabel( '$ \Delta \ln(%s^{\, sisL}/%s^{\, sisS})$' %(cell_vars[h], cell_vars[h]), fontsize= 21.0)
    if h == 1:
        ax.set_xlabel( '$|\\alpha_b| \, (V/\mu_V)$', fontsize= 21.0)
    else:
        ax.set_xlabel( '$|\\alpha_b| \, (%s)$' %(cell_vars[h]), fontsize= 23.0)
    #ax.text(0.48, 0.72, leg, fontsize=6, color = 'k', transform=ax.transAxes)
    ax.set_xlim([-0.1, 0.65])
    ax.set_ylim([-0.7, 0.4])
    ax.set_xticks([0.0, 0.3, 0.6])
    ax.set_yticks([-0.6, -0.3, 0.0, 0.3])
    fig.set_size_inches(5.0, 5.0)
    plt.gcf().subplots_adjust(bottom=0.3, left = 0.3)
    plt.savefig(resultsPath + "%s_alphab_Vs_log(SLdivSs)_%s.pdf"%(cell_vars[h],radius), format='pdf', dpi=200, transparent =True)
    plt.close(fig)


    print "NEW G T DATA ....."
    print len(all_avRelGrowthRatesCCTimes[0][h]), len(all_avRelGrowthRatesCCTimes[0][h] + all_avRelGrowthRatesCCTimes[1][h])
    print len(all_cellCycleTimes[0]), len(all_cellCycleTimes[0] + all_cellCycleTimes[1])

    print len(all_vars_birth[0][h])
    print len(all_vars_birth[0][h] + all_vars_birth[1][h])

    gTData = np.array(all_avRelGrowthRatesCCTimes[0][h] + all_avRelGrowthRatesCCTimes[1][h])
    allVarsBirth = np.array(all_vars_birth[0][h] + all_vars_birth[1][h]) 


    Asplit = sp_cell_size_birth
    Abins = [[] for k in range(len(Asplit)-1)]; dABins = [[] for k in range(len(Asplit)-1)]
    for k in range(len(allVarsBirth)):
        for l in range(len(Asplit)-1):
            if (allVarsBirth[k] >= Asplit[l] and allVarsBirth[k] < Asplit[l + 1]): 
                Abins[l].append(allVarsBirth[k]); dABins[l].append(gTData[k])   
    x = [np.median(Abins[k]) for k in range(len(Abins))]; y = [np.median(dABins[k]) for k in range(len(dABins))]
    yerr = [[y[k] - np.percentile(dABins[k], 25) for k in range(len(dABins))], [np.percentile(dABins[k], 75) - y[k] for k in range(len(dABins))]]           
    m, b = np.polyfit(x, y, 1); r, p = pearsonr(np.array(allVarsBirth), np.array(gTData)) 
    leg = "n = " + str(len(allVarsBirth)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b, 2)) + "\n"   + "c.c., p = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"    


    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(allVarsBirth, gTData, color= "c",  edgecolors='None', s=65, alpha = 0.1) 
    ax.plot(x, m * np.array(x) + b, linestyle = '-', color = 'r', linewidth = 1.0 )
    ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
    ax.set_ylabel( 'norm. $g_{rel} \\times T (%s)$' %(cell_vars[h]), fontsize= 21.0)
    if h == 1:
        ax.set_xlabel( 'birth size $(V/\mu_V)$', fontsize= 21.0)
    else:
        ax.set_xlabel( 'birth size $(%s)$' %(cell_vars[h]), fontsize= 23.0)
    ax.text(0.48, 0.72, leg, fontsize=7, color = 'k', transform=ax.transAxes)
    #ax.set_xlim([-0.6, 0.6])
    #ax.set_ylim([0.0, 1.6])
    #ax.set_xticks([-0.5, 0.0, 0.5])
    #ax.set_yticks([0.5, 1.0, 1.5])
    fig.set_size_inches(4.0, 4.0)
    plt.gcf().subplots_adjust(bottom=0.3, left = 0.3)
    plt.savefig(resultsPath + "%s_birthSize_Vs_gRelT_%s.pdf"%(cell_vars[h],radius), format='pdf', dpi=200, transparent =True)
    plt.close(fig)



    growthdata = np.array(all_av_rel_growth_sisters[0][h] + all_av_rel_growth_sisters[1][h])/np.mean(all_av_rel_growth_sisters[0][h] + all_av_rel_growth_sisters[1][h])
    asymdata = all_asym_sisters[0][h] + all_asym_sisters[1][h]

    Asplit = sp_cell_asym_birth
    Abins = [[] for k in range(len(Asplit)-1)]; dABins = [[] for k in range(len(Asplit)-1)]
    for k in range(len(growthdata)):
        for l in range(len(Asplit)-1):
            if (asymdata[k] >= Asplit[l] and asymdata[k] < Asplit[l + 1]): 
                Abins[l].append(asymdata[k]); dABins[l].append(growthdata[k])   
    x = [np.median(Abins[k]) for k in range(len(Abins))]; y = [np.median(dABins[k]) for k in range(len(dABins))]
    yerr = [[y[k] - np.percentile(dABins[k], 25) for k in range(len(dABins))], [np.percentile(dABins[k], 75) - y[k] for k in range(len(dABins))]]           
    m, b = np.polyfit(x, y, 1); r, p = pearsonr(np.array(asymdata), np.array(growthdata)) 
    leg = "n = " + str(len(asymdata)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b, 2)) + "\n"   + "c.c., p = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"    


    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(asymdata, growthdata, color= "c",  edgecolors='None', s=65, alpha = 0.1) 
    ax.plot(x, m * np.array(x) + b, linestyle = '-', color = 'r', linewidth = 1.0 )
    ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
    ax.set_ylabel( '$g_{rel}/$mean$(g_{rel}) (%s)$' %(cell_vars[h]), fontsize= 21.0)
    if h == 1:
        ax.set_xlabel( '$\\alpha_b \, (V/\mu_V)$', fontsize= 21.0)
    else:
        ax.set_xlabel( '$\\alpha_b \, (%s)$' %(cell_vars[h]), fontsize= 23.0)
    #ax.text(0.48, 0.72, leg, fontsize=6, color = 'k', transform=ax.transAxes)
    ax.set_xlim([-0.6, 0.6])
    ax.set_ylim([0.4, 1.6])
    ax.set_xticks([-0.5, 0.0, 0.5])
    ax.set_yticks([0.5, 1.0, 1.5])
    fig.set_size_inches(5.0, 5.0)
    plt.gcf().subplots_adjust(bottom=0.3, left = 0.3)
    plt.savefig(resultsPath + "%s_alphab_Vs_gRel_%s.pdf"%(cell_vars[h],radius), format='pdf', dpi=200, transparent =True)
    plt.close(fig)

    if h == 0: np.savetxt(resultsPath + "alphab_vs_gRel_vol_r_%s"%(radius) + ".out", [asymdata, growthdata], delimiter=",", fmt="%s")
    if h == 1: np.savetxt(resultsPath + "alphab_vs_gRel_normVol_r_%s"%(radius) + ".out", [asymdata, growthdata], delimiter=",", fmt="%s")

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(all_asym_birth, all_log_ratios_slopes, color= "c",  edgecolors='None', s=65, alpha = 0.1) 
    ax.plot(x, m * np.array(x) + b, linestyle = '-', color = 'r', linewidth = 1.0 )
    ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
    ax.set_ylabel( '$ \Delta \ln(%s^{\, sisL}/%s^{\, sisS})$' %(cell_vars[h], cell_vars[h]), fontsize= 21.0)
    if h == 1:
        ax.set_xlabel( '$|\\alpha_b| \, (V/\mu_V)$', fontsize= 21.0)
    else:
        ax.set_xlabel( '$|\\alpha_b| \, (%s)$' %(cell_vars[h]), fontsize= 23.0)
    #ax.text(0.48, 0.72, leg, fontsize=6, color = 'k', transform=ax.transAxes)
    ax.set_xlim([-0.1, 0.65])
    ax.set_ylim([-0.7, 0.4])
    ax.set_xticks([0.0, 0.3, 0.6])
    ax.set_yticks([-0.6, -0.3, 0.0, 0.3])
    fig.set_size_inches(4.0, 4.0)
    plt.gcf().subplots_adjust(bottom=0.3, left = 0.3)
    plt.savefig(resultsPath + "%s_alphab_Vs_log(SLdivSs)_%s.pdf"%(cell_vars[h],radius), format='pdf', dpi=200, transparent =True)
    plt.close(fig)



    Asplit = [ np.amin(all_ratios_birth), 1.3, 1.5, 1.8, np.amax(all_ratios_birth)]

    Abins = [[] for k in range(len(Asplit)-1)]; dABins = [[] for k in range(len(Asplit)-1)]
    for k in range(len(all_ratios_birth)):
        for l in range(len(Asplit)-1):
            if (all_ratios_birth[k] >= Asplit[l] and all_ratios_birth[k] < Asplit[l + 1]): 
                Abins[l].append(all_ratios_birth[k]); dABins[l].append(all_slopes[k])   
    x = [np.median(Abins[k]) for k in range(len(Abins))]; y = [np.median(dABins[k]) for k in range(len(dABins))]
    yerr = [[y[k] - np.percentile(dABins[k], 25) for k in range(len(dABins))], [np.percentile(dABins[k], 75) - y[k] for k in range(len(dABins))]]           
    m, b = np.polyfit(x, y, 1); r, p = pearsonr(np.array(all_ratios_birth), np.array(all_slopes)) 
    leg = "n = " + str(len(all_ratios_birth)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b, 2)) + "\n"   + "c.c., p = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"    

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    #ax.scatter(all_ratios_birth, all_slopes, c= newcolourCode,  edgecolors='None', alpha = 0.3) 
    ax.scatter(all_ratios_birth, np.array(all_slopes)/np.array(all_ratios_birth), color= "0.55",  edgecolors='None', alpha = 0.3) 
    ax.errorbar(x, y, yerr =  yerr, linestyle='-', color = 'k', marker='o')
    ax.set_ylabel( '$ d (%s^{\, sisL}/%s^{\, sisS})/ d (t / (T^{\, sisL} + T^{\, sisS}))$' %(cell_vars[h], cell_vars[h]), fontsize= 12.0)
    ax.set_xlabel( '${%s^{\, sisL}}_b / {%s^{\, sisS}}_b$' %(cell_vars[h], cell_vars[h]), fontsize= 10.0)
    ax.text(0.48, 0.72, leg, fontsize=6, color = 'k', transform=ax.transAxes)
    ax.set_xlim([0.75, 2.75])
    ax.set_ylim([-1.0, 0.5])
    ax.set_xticks([1.0, 1.5, 2.0, 2.5])
    ax.set_yticks([-0.8, -0.4, 0.0, 0.4])
    fig.set_size_inches(3.0, 2.2)
    plt.gcf().subplots_adjust(bottom=0.22, left = 0.2)
    plt.savefig(resultsPath + "%s_del_VLdivVs_sisters_%s.pdf"%(cell_vars[h],radius), format='pdf', dpi=200)
    plt.close(fig)


 

    all_ratios_birth = [all_ratio_sisters[h][i][0] - 1 for i in range(len(all_ratio_sisters[h])) if all_times_sisters[i][-1] != 0  ]
    all_slopes = [(all_ratio_sisters[h][i][-1] - all_ratio_sisters[h][i][0])/(all_times_sisters[i][-1] * all_ratio_sisters[h][i][0]) for i in range(len(all_ratio_sisters[h])) if all_times_sisters[i][-1] != 0 ] 

    Asplit = [ np.amin(all_ratios_birth), np.percentile(all_ratios_birth, 25), np.percentile(all_ratios_birth, 50), np.percentile(all_ratios_birth, 75), np.amax(all_ratios_birth)]
    Abins = [[] for k in range(len(Asplit)-1)]; dABins = [[] for k in range(len(Asplit)-1)]
    for k in range(len(all_ratios_birth)):
        for l in range(len(Asplit)-1):
            if (all_ratios_birth[k] >= Asplit[l] and all_ratios_birth[k] < Asplit[l + 1]): 
                Abins[l].append(all_ratios_birth[k]); dABins[l].append(all_slopes[k])   
    x = [np.median(Abins[k]) for k in range(len(Abins))]; y = [np.median(dABins[k]) for k in range(len(dABins))]
    yerr = [[y[k] - np.percentile(dABins[k], 25) for k in range(len(dABins))], [np.percentile(dABins[k], 75) - y[k] for k in range(len(dABins))]]           
    m, b = np.polyfit(x, y, 1); r, p = pearsonr(np.array(all_ratios_birth), np.array(all_slopes)) 
    leg = "n = " + str(len(all_ratios_birth)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "c.c., p = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"    

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(all_ratios_birth, all_slopes, c= newcolourCode,  edgecolors='None', alpha = 0.3) 
    ax.errorbar(x, y, yerr =  yerr, linestyle='-', color = 'k', marker='o')
    ax.set_ylabel( ' d (%s1 / %s2) / dt / (%s1 / %s2)' %(cell_vars[h], cell_vars[h], cell_vars[h], cell_vars[h]), fontsize= 8.0)
    ax.set_xlabel( '%s1 / %s2 - 1 at birth' %(cell_vars[h], cell_vars[h]), fontsize= 8.0)
    ax.text(0.07, 0.85, leg, fontsize=6, color = 'k', transform=ax.transAxes)
    ax.set_xlim([-0.5, 3.2])
    ax.set_ylim([-1.0, 0.5])
    plt.savefig(resultsPath + "%s_converg2_rates_sisters_%s.pdf"%(cell_vars[h],radius), format='pdf', dpi=200)
    plt.close(fig) 

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter((all_alpha_b[0][h] + all_alpha_b[1][h]), (all_alpha_ns_neigh[0] + all_alpha_ns_neigh[1]), color='c',  s = 80, marker ="o", edgecolors='None', alpha = 0.08) 
    m, b = np.polyfit((all_alpha_b[0][h] + all_alpha_b[1][h]), (all_alpha_ns_neigh[0] + all_alpha_ns_neigh[1]), 1); r, p = pearsonr((all_alpha_b[0][h] + all_alpha_b[1][h]), (all_alpha_ns_neigh[0] + all_alpha_ns_neigh[1])) 
    leg = "$N$ = " + str(len((all_alpha_b[0][h] + all_alpha_b[1][h]))) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
    ax.text(0.07, 0.73, leg, fontsize=13, color = 'k', transform=ax.transAxes)
    ax.axhline(y= 0.11, c='k')
    ax.axhline(y= -0.11, c='k')
    ax.axvline(x= 0.11, c='k')
    ax.axvline(x= -0.11, c='k')
    ax.set_ylabel( '$(V - {V^{ns-neigh}})/(V + {V^{ns-neigh}})$ ', fontsize= 19.0)
    ax.set_xlabel( '$\\alpha_b$', fontsize= 27.0)
    plt.gcf().subplots_adjust(bottom = 0.16, left = 0.16)
    fig.set_size_inches(6.5, 4.5)
    plt.savefig(resultsPath + "Asymm%sBirth_vs_AMeanNonSisNeighVBirth_%s.pdf"%(cell_vars[h], radius), format='pdf', dpi=200)
    plt.close(fig)


print "num_large_asym = ", num_large_asym, " num_small_asym = ", num_small_asym

print "birth_size_outlier_number ",  birth_size_outlier_number

print "birth_size_mid_number ",  birth_size_mid_number

