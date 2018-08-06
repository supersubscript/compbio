import cPickle
import numpy as np
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from scipy.stats.stats import skew
import all_plant_variables as apv 
from matplotlib.patches import Ellipse

# # # # # # # # # # # # # CELLS TO INCLUDE IN STATS: CENTRAL ZONE ONLY OR CLEANED FULL LINEAGES 

radius = "30"
zone = "all_data" # split data according to: "all_data", "early_times", "late_times", "inner_zone", "outer_zone", "small_volumes", "large_volumes", "random", "born_in_light", "born_in_dark", etc
resultsPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_pooled_sister_stats/" + zone + "/"
excludeSTD_div, excludeSTD_inc, excludeSTD_T, excludeSTD_sumsis = 20, 20, 20, 20

# # # # # # # # # # # # # VARIABLE LABELS 
cell_vars = ["V", "V_n", "A", "Aop",  "Aip", "Aa"]
cell_vars_tex = ["V", "{V/\mu_V}", "A", "A_{op}",  "A_{ip}", "A_a"]

# # # # # # # # # # # # # FIGURES AND FIGURE VARIABLES
plantFontSize = 7
plt.locator_params(axis = 'x', nbins=4); plt.locator_params(axis = 'y', nbins=5)
plt.rcParams['xtick.labelsize'] = 17; plt.rcParams['ytick.labelsize'] = 17

# # # # # # # # # # # # # FORMAT FOR EXPORTED DATA
ABirthVsADiv_fits = [[" var Bir. vs var Div.", "slope", "intercept", "corr. coeff.", "p-value", " var Bir. vs Inc.", "slope", "intercept", "corr. coeff.", "p-value", "mean var birth", "coeff. var. var birth", "mean var div.", "std var division/mean var birth", "std (var S1 - var S2)/ 2(var S1 + var S2)"]]

all_plants, all_cellCycleTimes, all_t_birth, all_distances, all_times, all_time_of_day_at_birth, all_time_of_day_at_division, all_vol_of_neighs_birth, all_asym_L1_vol_of_nonsis_neighs_birth = [], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]]
all_vars_birth, all_vars_div, all_vars_asym_div = [[[] for var in cell_vars], [[] for var in cell_vars]], [[[] for var in cell_vars], [[] for var in cell_vars]], [[[] for var in cell_vars], [[] for var in cell_vars]], 

# # # # # # # # # # # # # FUNCTION THAT SPLITS DATA
def include_sample(X):
    if(zone=="inner_zone"):
        if(dist_centre[0][X] + dist_centre[1][X] < med_dist_centre): return True
        else: return False
    elif(zone=="outer_zone"):
        if(dist_centre[0][X] + dist_centre[1][X] >= med_dist_centre): return True
        else: return False
    elif(zone=="early_times"):
        if(time_at_birth[0][X] < np.median(time_at_birth[0])): return True
        else: return False
    elif(zone=="first_quarter_times"):
        if(time_at_birth[0][X] < np.percentile(time_at_birth[0], 25)): return True
        else: return False
    elif(zone=="late_times"):
        if(time_at_birth[0][X] >= np.median(time_at_birth[0])): return True
        else: return False
    elif(zone=="small_mother_volumes"):
        if(size_vars_birth[0][0][X] + size_vars_birth[1][0][X]  < med_mother_vol): return True
        else: return False
    elif(zone=="large_mother_volumes"):
        if(size_vars_birth[0][0][X] + size_vars_birth[1][0][X]  >= med_mother_vol): return True
        else: return False
    elif(zone=="small_neighbour_vols"):
        if(vol_neighs_at_birth[0][X]  < np.median(vol_neighs_at_birth[0])): return True
        else: return False
    elif(zone=="large_neighbour_vols"):
        if(vol_neighs_at_birth[0][X]  >= np.median(vol_neighs_at_birth[0])): return True
        else: return False
    elif(zone=="small_L1_neighbour_vols"):
        if(L1_vol_neighs_at_birth[0][X]  < np.median(L1_vol_neighs_at_birth[0])): return True
        else: return False
    elif(zone=="large_L1_neighbour_vols"):
        if(L1_vol_neighs_at_birth[0][X]  >= np.median(L1_vol_neighs_at_birth[0])): return True
        else: return False
    elif(zone=="random"):
        if(np.random.rand() < 0.5): return True
        else: return False       
    elif(zone=="born_in_light"): #lights on between 0 and 16 hrs; lights off between 16 and 24hrs.
        if(time_of_day_at_birth[0][X] < 16): return True
        else: return False
    elif(zone=="born_in_dark"): 
        if(time_of_day_at_birth[0][X] >= 16 and time_of_day_at_birth[0][X] < 24 ): return True
        else: return False   
    elif(zone=="born_in_morning"): #lights on between 0 and 16 hrs; lights off between 16 and 24hrs.
        if(time_of_day_at_birth[0][X] < 8): return True
        else: return False   
    elif(zone=="born_in_afternoon"): #lights on between 0 and 16 hrs; lights off between 16 and 24hrs.
        if(time_of_day_at_birth[0][X] >= 8 and time_of_day_at_birth[0][X] < 16): return True
        else: return False  
    elif(zone=="born_in_late_afternoon"): #lights on between 0 and 16 hrs; lights off between 16 and 24hrs.
        if(time_of_day_at_birth[0][X] >= 12 and time_of_day_at_birth[0][X] < 16): return True
        else: return False
    elif(zone=="born_in_early_afternoon"): #lights on between 0 and 16 hrs; lights off between 16 and 24hrs.
        if(time_of_day_at_birth[0][X] >= 8 and time_of_day_at_birth[0][X] < 12): return True
    elif(zone=="born_around_midday"): 
        if(time_of_day_at_birth[0][X] >= 4 and time_of_day_at_birth[0][X] < 12): return True
        else: return False    
    elif(zone=="symmetric_birth"):
        if(size_vars_birth[0][0][i]/size_vars_birth[1][0][i] < 1.25 and size_vars_birth[0][0][i]/size_vars_birth[1][0][i] > 1.0/1.25): return True
        else: return False
    elif(zone=="asymmetric_birth"):
        if(size_vars_birth[0][0][i]/size_vars_birth[1][0][i] > 1.25): return True
        elif(size_vars_birth[0][0][i]/size_vars_birth[1][0][i] < 1.0/1.25): return True
        else: return False
    else: return True

def round_sig(x, sig=2):
    return np.round(x, sig-int(np.floor(np.log10(x)))-1)

def outlier(x, y, xp, yp , std):
    lambda_, v = np.linalg.eig(np.cov(x,y)); lambda_ = np.sqrt(lambda_) 
    theta = -np.arccos(v[0, 0])
    x0 = xp - np.mean(x); y0 = yp- np.mean(y)
    x1 =  np.cos(theta)*x0 - np.sin(theta)*y0; y1 =  np.sin(theta)*x0 + np.cos(theta)*y0
    if x1/lambda_[0] * x1/lambda_[0] + y1/lambda_[1] * y1/lambda_[1] < std * std: return False
    else : return True
    

# # # # # # # # # # # # # MAIN LOOP
for plant in range(len(apv.allPlants)):
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
    time_at_birth, time_of_day_at_birth, time_of_day_at_division, vol_neighs_at_birth, L1_vol_neighs_at_birth, L1_mean_vol_nonsis_neighs_at_birth  = [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]]    
    cell_cycle_time = [[],[]]

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
                d0 = data_t0[sisInd[k]][2]
                distance_over_cycle[k].append([]); time_since_birth[k].append([]); dist_evolution = []                
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
                        time_of_day_at_division[k].append(times[data_t1[sisInd[k]][0] + 1]%24)
                        neigh1 = dataSeg[0]['neigbourhood'][data_t0[sisInd[0]][2]]; neigh2 = dataSeg[0]['neigbourhood'][data_t0[sisInd[1]][2]]; exc = [x for x in neigh1 if x not in neigh2]                                 
                        all_neighs = [x for x in exc + neigh2 if x not in [data_t0[sisInd[0]][2], data_t0[sisInd[1]][2]]]  
                        vol_neighs_at_birth[k].append(sum([dataSeg[0]['volumes'][x] for x in all_neighs if x != 1])/meanVolumes[j])
                        L1_neighs = [x for x in exc+neigh2 if 1 in dataSeg[0]['neigbourhood'][x] ]
                        L1_vol_neighs_at_birth[k].append(sum([dataSeg[0]['volumes'][x] for x in L1_neighs])/meanVolumes[j])
                        non_sis_neighs = [x for x in dataSeg[0]['neigbourhood'][d0] if x != sisInd[(k + 1)%2]] 
                        non_sis_L1_neighs = [x for x in non_sis_neighs if 1 in dataSeg[0]['neigbourhood'][x] ]
                        L1_mean_vol_nonsis_neighs_at_birth[k].append(sum([dataSeg[0]['volumes'][x] for x in non_sis_L1_neighs])/len(non_sis_L1_neighs))

                    size_vars[k][0][1][-1].append(dataSeg[0]['volumes'][d0]); size_vars[k][1][1][-1].append(dataSeg[0]['volumes'][d0]/meanVolumes[j])  
                    time_since_birth[k][-1].append(times[j] - times[data_t0[sisInd[k]][0] + 1])

                    if(time_since_birth[k][-1][-1]< 0): 
                        print "BUG  ... ", j, data_t0[sisInd[k]], time_since_birth[-1] 

                    xSq = (dataSeg[0]['barycenter'][d0][0] - L1centre[j][0]) * (dataSeg[0]['barycenter'][d0][0] - L1centre[j][0])
                    ySq = (dataSeg[0]['barycenter'][d0][1] - L1centre[j][1]) * (dataSeg[0]['barycenter'][d0][1] - L1centre[j][1])
                    zSq = (dataSeg[0]['barycenter'][d0][2] - L1centre[j][2]) * (dataSeg[0]['barycenter'][d0][2] - L1centre[j][2])
                    dist_evolution.append(np.sqrt(xSq + ySq + zSq))
                    distance_over_cycle[k][-1].append(np.sqrt(xSq + ySq + zSq)) 

                    pArea = 0; sArea = 0; oPArea = 0
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
                    newDaughters = [x[0] for x in dataLineage if x[1] == d0]
                    if (j != data_t1[sisInd[k]][0]) and (len(newDaughters) != 1) : print " More or less than one offspring ...  "
                    if len(newDaughters) == 0 : print " Zero offspring ...  "; break
                    else : d0 = newDaughters[0] 

                dist_centre[k].append(np.mean(dist_evolution))
                fobj = file(segDataPathList[data_t1[sisInd[k]][0] + 1])
                dataSeg = cPickle.load(fobj)
                
                vol, volNorm, sArea, pArea, oPArea = [0 for j in range(len(newDaughters))], [0 for j in range(len(newDaughters))], [0 for j in range(len(newDaughters))], [0 for j in range(len(newDaughters))], [0 for j in range(len(newDaughters))]

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

#                if len(newDaughters) == 2 : 
#                    random = 1 if np.random.rand() < 0.5 else -1
#                    asymmetry_div[k][0].append(0.5*(vol[0] - vol[1])*random/np.sum(vol)); asymmetry_div[k][1].append(0.5*(volNorm[0] - volNorm[1])*random/np.sum(volNorm));  asymmetry_div[k][2].append(0.5*(sArea[0] - sArea[1])*random/np.sum(sArea));  asymmetry_div[k][4].append(0.5* random*(oPArea[0] -oPArea[1])/np.sum(oPArea)); asymmetry_div[k][5].append(0.5* random*(sArea[0] - sArea[1] - (pArea[0] - pArea[1]) -  (oPArea[0] - oPArea[1]))/(np.sum(sArea) - np.sum(pArea) - np.sum(oPArea))) 
#                    if np.sum(pArea) != 0 : asymmetry_div[k][3].append(0.5* random*(pArea[0] - pArea[1])/np.sum(pArea))
#                    else : asymmetry_div[k][3].append(0)
#                else: 
#                    for asymmetry in asymmetry_div[k]: asymmetry.append(0)
#                    print " Aberant number of daughters =  ", len(newDaughters)   


    size_vars_birth = [[],[]]
    for var in size_vars[0]: size_vars_birth[0].append(np.array([var[1][i][0] for i in range(0, len(var[1]))]))
    for var in size_vars[1]: size_vars_birth[1].append(np.array([var[1][i][0] for i in range(0, len(var[1]))]))

    print np.array(size_vars_birth[0][0])/np.array(size_vars_birth[1][0])
    print np.array(size_vars_birth[0][1])/np.array(size_vars_birth[1][1])

    med_mother_vol = np.median(np.array(size_vars_birth[0][0]) + np.array(size_vars_birth[1][0]))
    med_dist_centre = np.median(np.array(dist_centre[0]) + np.array(dist_centre[1]))

    to_include = []
    for i in range(0, len(cell_cycle_time[0])):      
         if(include_sample(i)): to_include.append(i)
    print "to include: ", len(to_include), len(cell_cycle_time[0])

    mean_vars_birth = []

    for var_len in range(len(size_vars_birth[0])): mean_vars_birth.append(np.mean([size_vars_birth[0][var_len][j] for j in to_include] + [size_vars_birth[1][var_len][j] for j in to_include])) 
    mean_vars_birth[1] = 1.0;  #normalized cell volume stats
    mean_vars_div = []
    for var_len in range(len(size_vars_div[0])): mean_vars_div.append(np.mean([size_vars_div[0][var_len][j] for j in to_include] + [size_vars_div[1][var_len][j] for j in to_include]))
    mean_period = np.median([cell_cycle_time[0][j] for j in to_include] + [cell_cycle_time[1][j] for j in to_include])  

    for i in to_include:      
        all_plants.append(plant); 
        for k in range(2):
            all_time_of_day_at_division[k].append(time_of_day_at_division[k][i])
            all_time_of_day_at_birth[k].append(time_of_day_at_birth[k][i])
            all_vol_of_neighs_birth[k].append(vol_neighs_at_birth[k][i])
            all_asym_L1_vol_of_nonsis_neighs_birth[k].append((size_vars_birth[k][0][i] - L1_mean_vol_nonsis_neighs_at_birth[k][i])/(L1_mean_vol_nonsis_neighs_at_birth[k][i]+size_vars_birth[k][0][i]))
            all_distances[k].append(dist_tbirth[k][i]/np.power(mean_vars_birth[0], 0.3333333))
            all_cellCycleTimes[k].append(cell_cycle_time[k][i]/mean_period)
            all_t_birth[k].append(time_at_birth[k][i])           

            for j in range(len(cell_vars)):                
                all_vars_birth[k][j].append(size_vars_birth[k][j][i]/mean_vars_birth[j])  
                all_vars_div[k][j].append(size_vars_div[k][j][i]/mean_vars_birth[j])  

    print np.array(all_vars_birth[0][0])/np.array(all_vars_birth[1][0]) 
    print np.array(all_vars_birth[0][1])/np.array(all_vars_birth[1][1])  


print len(all_vars_birth[0][0]), len(all_vars_birth[1][0]) 
print len(all_vars_div[0][0]), len(all_vars_div[1][0]) 
print len(all_cellCycleTimes[0]), len(all_cellCycleTimes[1]) 
print len(all_vol_of_neighs_birth[0]), len(all_vol_of_neighs_birth[1]) 


ratio_distances = np.array(all_distances[0])/np.array(all_distances[1])
ratio_birth_vols = np.array(all_vars_birth[0][0] )/np.array(all_vars_birth[1][0])
ratios_div_vol_birth_vol = (np.array(all_vars_div[0][0]) + np.array(all_vars_div[1][0]))/(np.array(all_vars_birth[0][0]) + np.array(all_vars_birth[0][1]))


plt.hist(ratios_div_vol_birth_vol, bins =6)
plt.xlabel('$V_d^{\, sis}/V_b^{\, sis}$')
plt.savefig(resultsPath + "hist_ratio_div_birth_vol_%s.pdf"%radius, format='pdf', dpi=300)
plt.close()


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.scatter(ratio_distances,ratio_birth_vols, c= all_plants, marker ="o", edgecolors='None', alpha = 0.3) 
m, b = np.polyfit(ratio_distances,ratio_birth_vols, 1); r, p = pearsonr(ratio_distances,ratio_birth_vols) 
ax.plot(ratio_distances, m * ratio_distances + b, linestyle = '-', color = 'b' )
leg = "n = " + str(len(ratio_distances)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "c.c., p = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
ax.text(0.07, 0.85, leg, fontsize=6, color = 'k', transform=ax.transAxes)
ax.set_xlabel( 'ratio sister distances from centre ', fontsize= 8.0)
ax.set_ylabel( 'ratio sister volumes ', fontsize= 8.0)
plt.savefig(resultsPath + "sisDistRatio_Vs_VolRatio_r_%s.pdf"%(radius), format='pdf', dpi=200)
plt.close(fig)


for k in range(len(cell_vars)):
    S_birth, S_div, Inc, T_s = [[], []], [[], []], [[], []], [[], []]
    A_birth, A_div, A_inc, A_T = [], [], [], []
    sis_dist, sis_time_birth = [], []
    Sum_sis_birth, Sum_sis_div = [], [] 
    A_non_sis_neigh_birth = [[], []]

    for i in range(len(all_vars_birth[0][0])):
        random = np.random.rand() 
    
        if(random < 0.5):  s1 = 0; s2 = 1

        else: s1 = 1; s2 = 0

        if((all_vars_birth[s1][k][i] + all_vars_birth[s2][k][i]) != 0 and (all_vars_div[s1][k][i] + all_vars_div[s2][k][i]) != 0 and  (all_vars_div[s1][k][i] + all_vars_div[s2][k][i] - all_vars_birth[s1][k][i] - all_vars_birth[s2][k][i]) != 0 and (all_cellCycleTimes[s1][i] + all_cellCycleTimes[s2][i]) != 0) :
            #print all_vars_birth[s1][k][i]/all_vars_birth[s2][k][i]
            A_birth.append((all_vars_birth[s1][k][i] - all_vars_birth[s2][k][i])/(all_vars_birth[s1][k][i] + all_vars_birth[s2][k][i]))
            A_div.append((all_vars_div[s1][k][i] - all_vars_div[s2][k][i])/(all_vars_div[s1][k][i] + all_vars_div[s2][k][i]))
            #A_div.append((all_vars_div[s1][k][i] - all_vars_div[s2][k][i])/(all_vars_birth[s1][k][i] + all_vars_birth[s2][k][i]))
            A_inc.append((all_vars_div[s1][k][i] - all_vars_birth[s1][k][i] + all_vars_birth[s2][k][i] - all_vars_div[s2][k][i])/(all_vars_div[s1][k][i] + all_vars_div[s2][k][i] - all_vars_birth[s1][k][i] - all_vars_birth[s2][k][i]))
            #A_inc.append((all_vars_div[s1][k][i] - all_vars_birth[s1][k][i] + all_vars_birth[s2][k][i] - all_vars_div[s2][k][i])/(all_vars_birth[s1][k][i] + all_vars_birth[s2][k][i]))

            A_T.append((all_cellCycleTimes[s1][i] - all_cellCycleTimes[s2][i])/(all_cellCycleTimes[s1][i] + all_cellCycleTimes[s2][i]))

            S_birth[0].append(all_vars_birth[s1][k][i]); S_birth[1].append(all_vars_birth[s2][k][i])
            S_div[0].append(all_vars_div[s1][k][i]); S_div[1].append(all_vars_div[s2][k][i])
            Inc[0].append(all_vars_div[s1][k][i] - all_vars_birth[s1][k][i]); Inc[1].append(all_vars_div[s2][k][i] - all_vars_birth[s2][k][i])
            T_s[0].append(all_cellCycleTimes[s1][i]); T_s[1].append(all_cellCycleTimes[s2][i])   

            sis_dist.append((all_distances[s1][i] + all_distances[s2][i])*0.5)  
            sis_time_birth.append(all_t_birth[s1][i])  
            Sum_sis_birth.append((all_vars_birth[s1][k][i] + all_vars_birth[s2][k][i]))
            Sum_sis_div.append((all_vars_div[s1][k][i] + all_vars_div[s2][k][i]))

            A_non_sis_neigh_birth[0].append(all_asym_L1_vol_of_nonsis_neighs_birth[s1][i]) ; A_non_sis_neigh_birth[1].append(all_asym_L1_vol_of_nonsis_neighs_birth[s2][i])

    print cell_vars[k], " median, mean ratio (S_b + S_b^sis)/ (S_d + S_d^sis) "
    print len(S_div[0]), len(np.array(S_birth[0]) + np.array(S_birth[1]))
    print np.median((np.array(S_birth[0]) + np.array(S_birth[1]))/(np.array(S_div[0]) + np.array(S_div[1]))), np.mean((np.array(S_birth[0]) + np.array(S_birth[1]))/(np.array(S_div[0]) + np.array(S_div[1])))

    print cell_vars[k], " median, mean ratio (S_b + S_b^sis)/ (delta + delta^sis) "
    print np.median((np.array(S_birth[0]) + np.array(S_birth[1]))/(np.array(S_div[0]) + np.array(S_div[1]) - np.array(S_birth[0]) - np.array(S_birth[1]))), np.mean((np.array(S_birth[0]) + np.array(S_birth[1]))/(np.array(S_div[0]) + np.array(S_div[1]) - np.array(S_birth[0]) - np.array(S_birth[1])))
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
   # ax.scatter(np.array(S_birth[0]) + np.array(S_birth[1]), np.absolute(np.array(A_birth)), c= all_plants, marker ="o", edgecolors='None', alpha = 0.3) 
    ax.scatter(np.array(S_birth[0]) + np.array(S_birth[1]), np.absolute(np.array(A_birth)), color="0.3", marker ="o", edgecolors='None', alpha = 0.3) 
    m, b = np.polyfit(np.array(S_birth[0]) + np.array(S_birth[1]), np.absolute(np.array(A_birth)), 1); r, p = pearsonr(np.array(S_birth[0]) + np.array(S_birth[1]), np.absolute(np.array(A_birth))) 
    #ax.plot(np.array(S_birth[0]) + np.array(S_birth[1]), m * (np.array(S_birth[0]) + np.array(S_birth[1])) + b, linestyle = '-', color = 'b' )
    leg = "n = " + str(len(S_div[0])) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "c.c., p = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
    ax.text(0.07, 0.85, leg, fontsize=8, color = 'k', transform=ax.transAxes)
    ax.set_xlabel( '$%s$ upon mother division'%cell_vars[k], fontsize= 14.0)
    ax.set_ylabel( '$|\\alpha_b|$', fontsize = 20.0)
    plt.gcf().subplots_adjust(bottom = 0.16, left = 0.13)
    fig.set_size_inches(6.0, 5.0)
    plt.savefig(resultsPath + "%sDiv_Vs_Asymm%s.pdf"%(cell_vars[k], radius), format='pdf', dpi=200)
   
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(np.array(A_birth)*np.array(A_birth), (np.array(S_birth[0]) + np.array(S_birth[1]))/(np.array(S_div[0]) + np.array(S_div[1])), color='c', marker ="o", s =80, edgecolors='None', alpha = 0.15) 
    m, b = np.polyfit(np.array(A_birth)*np.array(A_birth), (np.array(S_birth[0]) + np.array(S_birth[1]))/(np.array(S_div[0]) + np.array(S_div[1])), 1); r, p = pearsonr(np.array(A_birth)*np.array(A_birth), (np.array(S_birth[0]) + np.array(S_birth[1]))/(np.array(S_div[0]) + np.array(S_div[1]))) 
    leg = "$N$ = " + str(len(S_div[0])) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
    print leg
    leg_two = "$(%s)$"%cell_vars_tex[k]
    ax.text(0.07, 0.75, leg, fontsize=13, color = 'k', transform=ax.transAxes)
    ax.text(0.75, 0.8, leg_two, fontsize=25, color = 'k', transform=ax.transAxes)
    if k == 1:
        ax.set_ylabel( 'sis. sum at birth/sis. sum at div.', fontsize= 18.0)
        #ax.set_ylabel( '$({(V/\mu_V)}_b + {{(V/\mu_V)}_b}^{\, sis})/({(V/\mu_V)}_d + {{(V/\mu_V)}_d}^{\, sis})$', fontsize= 18.0)
    else:
        ax.set_ylabel( 'sis. sum at birth/sis. sum at div.', fontsize= 18.0)
        #ax.set_ylabel( '$({%s}_b + {{%s}_b}^{\, sis})/({%s}_d + {{%s}_d}^{\, sis})$'%(cell_vars_tex[k],cell_vars_tex[k],cell_vars_tex[k],cell_vars_tex[k]), fontsize= 18.0)
    ax.set_xlabel( '${\\alpha_b}^2$', fontsize = 28.0)
    plt.gcf().subplots_adjust(bottom = 0.18, left = 0.16)
    fig.set_size_inches(7.0, 6.0)
    plt.savefig(resultsPath + "%s_SumSisBirthsDivSumSisDivs_Vs_AsymmSq_%s.pdf"%(cell_vars[k], radius), format='pdf', dpi=200)
    plt.close(fig)


    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(sis_dist, np.absolute(np.array(A_birth)), color= 'c',  s =60, alpha = 0.2,  marker ="o", edgecolors='None') 
    m, b = np.polyfit(sis_dist, np.absolute(np.array(A_birth)), 1); r, p = pearsonr(sis_dist, np.absolute(np.array(A_birth))) 
    leg = "$N$ = " + str(len(sis_dist)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
    varTex = "($%s$)"%cell_vars_tex[k]
    ax.text(0.75, 0.85, varTex, fontsize=30, color = 'k', transform=ax.transAxes)
    ax.text(0.07, 0.75, leg, fontsize=13, color = 'k', transform=ax.transAxes)
    ax.set_xlabel( 'norm. distance from O', fontsize= 20.0)
    ax.set_ylabel( '$|\\alpha_b|$', fontsize= 30.0)
    plt.gcf().subplots_adjust(bottom = 0.16, left = 0.19)
    fig.set_size_inches(6.0, 5.0)
    plt.savefig(resultsPath + "dist_Vs_%sBirth_Asymm_%s.pdf"%(cell_vars[k], radius), format='pdf', dpi=200)
    plt.close(fig)


    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(sis_dist, np.absolute(np.array(A_T)), color= 'c',  s =60, alpha = 0.2,  marker ="o", edgecolors='None') 
    m, b = np.polyfit(sis_dist, np.absolute(np.array(A_T)), 1); r, p = pearsonr(sis_dist, np.absolute(np.array(A_T))) 
    leg = "$N$ = " + str(len(sis_dist)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
    ax.text(0.07, 0.85, leg, fontsize=13, color = 'k', transform=ax.transAxes)
    ax.set_xlabel( 'norm. distance from $O$', fontsize= 20.0)
    ax.set_ylabel( '$|\\alpha_T|$', fontsize= 30.0)
    plt.gcf().subplots_adjust(bottom = 0.16, left = 0.19)
    fig.set_size_inches(6.0, 5.0)
    plt.savefig(resultsPath + "dist_Vs_T_Asymm_%s.pdf"%(radius), format='pdf', dpi=200)
    plt.close(fig)

    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(sis_dist, np.absolute(np.array(A_div)), color= 'c',  s =60, alpha = 0.2, marker ="o", edgecolors='None') 
    m, b = np.polyfit(sis_dist, np.absolute(np.array(A_div)), 1); r, p = pearsonr(sis_dist, np.absolute(np.array(A_div))) 
    varTex = "($%s$)"%cell_vars_tex[k]
    ax.text(0.75, 0.85, varTex, fontsize=30, color = 'k', transform=ax.transAxes)
    leg = "$N$ = " + str(len(sis_dist)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
    ax.text(0.07, 0.75, leg, fontsize=13, color = 'k', transform=ax.transAxes)
    ax.set_xlabel( 'norm. distance from O', fontsize= 20.0)
    ax.set_ylabel( '$|\\alpha_d|$', fontsize= 30.0)
    plt.gcf().subplots_adjust(bottom = 0.16, left = 0.19)
    fig.set_size_inches(6.0, 5.0)
    plt.savefig(resultsPath + "dist_Vs_%sDiv_Asymm_%s.pdf"%(cell_vars[k], radius), format='pdf', dpi=200)
    plt.close(fig)

    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(sis_time_birth, np.absolute(np.array(A_birth)), color= 'c',  s =60, alpha = 0.15, marker ="o", edgecolors='None') 
    m, b = np.polyfit(sis_time_birth, np.absolute(np.array(A_birth)), 1); r, p = pearsonr(sis_time_birth, np.absolute(np.array(A_birth))) 
    varTex = "($%s$)"%cell_vars_tex[k]
    ax.text(0.75, 0.85, varTex, fontsize=30, color = 'k', transform=ax.transAxes)
    leg = "$N$ = " + str(len(sis_time_birth)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
    ax.text(0.07, 0.75, leg, fontsize=13, color = 'k', transform=ax.transAxes)
    ax.set_xlabel( 'time (h)', fontsize= 20.0)
    ax.set_ylabel( '$|\\alpha_b|$', fontsize= 30.0)
    plt.gcf().subplots_adjust(bottom = 0.16, left = 0.19)
    fig.set_size_inches(6.0, 5.0)
    plt.savefig(resultsPath + "time_Vs_%sBirth_Asymm_%s.pdf"%(cell_vars[k], radius), format='pdf', dpi=200)
    plt.close(fig)


    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(sis_time_birth, np.absolute(np.array(A_T)),  color= 'c',  s =60, alpha = 0.2, marker ="o", edgecolors='None') 
    m, b = np.polyfit(sis_time_birth, np.absolute(np.array(A_T)), 1); r, p = pearsonr(sis_time_birth, np.absolute(np.array(A_T))) 
    leg = "$N$ = " + str(len(sis_time_birth)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
    ax.text(0.07, 0.75, leg, fontsize=13, color = 'k', transform=ax.transAxes)
    ax.set_xlabel( 'time (h)', fontsize= 20.0)
    ax.set_ylabel( '$|\\alpha_T|$', fontsize= 30.0)
    plt.gcf().subplots_adjust(bottom = 0.16, left = 0.19)
    fig.set_size_inches(6.0, 5.0)
    plt.savefig(resultsPath + "time_Vs_T_Asymm_%s.pdf"%(radius), format='pdf', dpi=200)
    plt.close(fig)

    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(sis_time_birth, np.absolute(np.array(A_div)), color= 'c',  s =60, alpha = 0.2, marker ="o", edgecolors='None') 
    m, b = np.polyfit(sis_time_birth, np.absolute(np.array(A_div)), 1); r, p = pearsonr(sis_time_birth, np.absolute(np.array(A_div))) 
    leg = "$N$ = " + str(len(sis_time_birth)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
    varTex = "($%s$)"%cell_vars_tex[k]
    ax.text(0.75, 0.85, varTex, fontsize=30, color = 'k', transform=ax.transAxes)
    ax.text(0.07, 0.75, leg, fontsize=13, color = 'k', transform=ax.transAxes)
    ax.set_xlabel( 'time (h)', fontsize= 20.0)
    ax.set_ylabel( '$|\\alpha_d|$', fontsize= 30.0)
    plt.gcf().subplots_adjust(bottom = 0.16, left = 0.19)
    fig.set_size_inches(6.0, 5.0)
    plt.savefig(resultsPath + "time_Vs_%sDiv_Asymm_%s.pdf"%(cell_vars[k], radius), format='pdf', dpi=200)
    plt.close(fig)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    AminB = [-x for x in A_birth]
    print len(S_birth[0] + S_birth[1])
    print len(A_birth + AminB)

    vec = S_birth[0] + S_birth[1]
    if k == 1: vec = np.array(vec)/np.mean(np.array(vec))
    ax.scatter(A_birth + AminB, vec, marker ="o", color= 'c',  s =80, alpha = 0.1, edgecolors='None') 
    m, b = np.polyfit(A_birth + AminB, vec, 1); r, p = pearsonr(A_birth + AminB, vec) 
    leg = "$N$ = " + str(len(vec)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
    ax.text(0.07, 0.73, leg, fontsize=13, color = 'k', transform=ax.transAxes)
    if k ==1 :
        ax.set_ylabel( 'norm. birth vol./mean(norm. birth vol.)', fontsize= 19.0)
    else:
        ax.set_ylabel( '${%s}_b/\mu_b$'%cell_vars[k], fontsize= 19.0)        
    ax.set_xlabel( '$\\alpha_b$', fontsize= 27.0)
    ax.axhline(y= 1.16, c='k')
    ax.axhline(y= 0.84, c='k')
    ax.axvline(x= 0.11, c='k')
    ax.axvline(x= -0.11, c='k')
    plt.gcf().subplots_adjust(bottom = 0.16, left = 0.16)
    fig.set_size_inches(6.5, 4.5)
    plt.savefig(resultsPath + "Asymm%sBirth_vs_%sBirth_%s.pdf"%(cell_vars[k], cell_vars[k], radius), format='pdf', dpi=200)
    plt.close(fig)

    if k == 1: np.savetxt(resultsPath + "alphab_vs_NormbirthVolOverMeanNormbirthVol_r_%s"%(radius) + ".out", [A_birth + AminB, vec], delimiter=",", fmt="%s")

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    vec = A_non_sis_neigh_birth[0] + A_non_sis_neigh_birth[1]
    print 'A_non_sis: ', len(vec)
    print 'A-birth-sis ', len(A_birth + AminB)
    ax.scatter(A_birth + AminB, vec,   marker ="o",color= 'c',  s =60, alpha = 0.1, edgecolors='None') 
    m, b = np.polyfit(A_birth + AminB, vec, 1); r, p = pearsonr(A_birth + AminB, vec) 
    leg = "n = " + str(len(vec)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "c.c., p = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
    ax.text(0.07, 0.78, leg, fontsize=8, color = 'k', transform=ax.transAxes)
    ax.set_ylabel( '$(V_b - {V^{ns-neigh}}_b)/(V_b + {V^{ns-neigh}}_b)$', fontsize= 14.0)
    ax.set_xlabel( '$\\alpha_b$', fontsize= 20.0)
    plt.gcf().subplots_adjust(bottom = 0.16, left = 0.16)
    fig.set_size_inches(5.0, 4.0)
    plt.savefig(resultsPath + "Asymm%sBirth_vs_AMeanNonSisNeighVBirth_%s.pdf"%(cell_vars[k], radius), format='pdf', dpi=200)
    plt.close(fig)


    Asplit = [ np.amin(S_div[0]), np.percentile(S_div[0], 33.3),  np.percentile(S_div[0], 66.7),  np.amax(S_div[0])]
    Abins = [[] for h in range(len(Asplit)-1)]
    dABins = [[] for h in range(len(Asplit)-1)]
    for h in range(len(S_div[0])):
        for l in range(len(Asplit)-1):
            if (S_div[0][h] >= Asplit[l] and S_div[0][h] < Asplit[l + 1]): 
                Abins[l].append(S_div[0][h]); dABins[l].append(S_div[1][h])   
    x = [np.median(Abins[h]) for h in range(len(Abins))]; y = [np.median(dABins[h]) for h in range(len(dABins))]
    yerr = [[y[h] - np.percentile(dABins[h], 25) for h in range(len(dABins))], [np.percentile(dABins[h], 75) - y[h] for h in range(len(dABins))]]    

    fig = plt.figure()
    ax = fig.add_subplot(2, 3, 1)
    ax.scatter(S_div[0], S_div[1], c= all_plants, marker ="o", edgecolors='None', color= 'c',  s =60, alpha = 0.1) 
#ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
    m, b = np.polyfit(S_div[0], S_div[1], 1); r, p = pearsonr(S_div[0], S_div[1]) 
    ax.plot(S_div[0], m * np.array(S_div[0]) + b, linestyle = '-', color = 'b' )
    ax.errorbar(x, y, yerr =  yerr, linestyle='-', color = 'k', marker='o')
    leg = "n = " + str(len(S_div[0])) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "c.c., p = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
    ax.text(0.07, 0.85, leg, fontsize=6, color = 'k', transform=ax.transAxes)
    ax.set_xlabel( '%s at division, sister 1'%cell_vars[k], fontsize= 8.0)
    ax.set_ylabel( '%s at division, sister 2'%cell_vars[k], fontsize= 8.0)

    Asplit = [ np.amin(Inc[0]), np.percentile(Inc[0], 33.3),  np.percentile(Inc[0], 66.7), np.amax(Inc[0])]
    Abins = [[] for h in range(len(Asplit)-1)]
    dABins = [[] for h in range(len(Asplit)-1)]
    for h in range(len(Inc[0])):
        for l in range(len(Asplit)-1):
            if (Inc[0][h] >= Asplit[l] and Inc[0][h] < Asplit[l + 1]): 
                Abins[l].append(Inc[0][h]); dABins[l].append(Inc[1][h])   
    x = [np.median(Abins[h]) for h in range(len(Abins))]; y = [np.median(dABins[h]) for h in range(len(dABins))]
    yerr = [[y[h] - np.percentile(dABins[h], 25) for h in range(len(dABins))], [np.percentile(dABins[h], 75) - y[h] for h in range(len(dABins))]] 
    ax = fig.add_subplot(2, 3, 2)
    ax.scatter(Inc[0], Inc[1], c= all_plants, marker ="o", edgecolors='None', color= 'c',  s =60, alpha = 0.1) 
    #ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
    m, b = np.polyfit(Inc[0], Inc[1], 1); r, p = pearsonr(Inc[0], Inc[1]) 
    ax.plot(Inc[0], m * np.array(Inc[0]) + b, linestyle = '-', color = 'b' )
    ax.errorbar(x, y, yerr =  yerr, linestyle='-', color = 'k', marker='o')
    leg = "n = " + str(len(Inc[0])) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "c.c., p = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
    ax.text(0.07, 0.85, leg, fontsize=6, color = 'k', transform=ax.transAxes)
    ax.set_xlabel( '%s increment, sister 1'%cell_vars[k], fontsize= 8.0)
    ax.set_ylabel( '%s increment, sister 2'%cell_vars[k], fontsize= 8.0)

    Asplit = [ np.amin(T_s[0]), np.percentile(T_s[0], 33.3),  np.percentile(T_s[0], 66.7), np.amax(T_s[0])]
    Abins = [[] for h in range(len(Asplit)-1)]
    dABins = [[] for h in range(len(Asplit)-1)]
    for h in range(len(Inc[0])):
        for l in range(len(Asplit)-1):
            if (T_s[0][h] >= Asplit[l] and T_s[0][h] < Asplit[l + 1]): 
                Abins[l].append(T_s[0][h]); dABins[l].append(T_s[1][h])   
    x = [np.median(Abins[h]) for h in range(len(Abins))]; y = [np.median(dABins[h]) for h in range(len(dABins))]
    #yerr = [[y[h] - np.percentile(dABins[h], 25) for h in range(len(dABins))], [np.percentile(dABins[h], 75) - y[h] for h in range(len(dABins))]] 
    ax = fig.add_subplot(2, 3, 3)
    ax.scatter(T_s[0], T_s[1], c= all_plants, marker ="o", edgecolors='None', color= 'c',  s =60, alpha = 0.1) 
    #ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
    m, b = np.polyfit(T_s[0], T_s[1], 1); r, p = pearsonr(T_s[0], T_s[1]) 
    ax.plot(T_s[0], m * np.array(T_s[0]) + b, linestyle = '-', color = 'b' )
    #ax.errorbar(x, y, yerr =  yerr, linestyle='-', color = 'k', marker='o')
    leg = "n = " + str(len(T_s[0])) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "c.c., p = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
    ax.text(0.07, 0.85, leg, fontsize=6, color = 'k', transform=ax.transAxes)
    ax.set_xlabel( 'Cell cycle time, sister 1 ', fontsize= 8.0)
    ax.set_ylabel( 'Cell cycle time, sister 2', fontsize= 8.0)
   
    A_birth_no = [A_birth[i] for i in range(len(A_birth)) if not outlier(A_birth, A_div, A_birth[i], A_div[i], excludeSTD_div)]
    A_div_no = [A_div[i] for i in range(len(A_div)) if not outlier(A_birth, A_div, A_birth[i], A_div[i], excludeSTD_div)]
    all_plants_no = [all_plants[i] for i in range(len(all_plants)) if not outlier(A_birth, A_div, A_birth[i], A_div[i], excludeSTD_div)]

    if k == 1: np.savetxt(resultsPath + "alphab_vs_alphad_r_%s"%(radius) + ".out", [A_birth_no, A_div_no], delimiter=",", fmt="%s")

    Asplit = [ np.amin(A_birth_no), np.percentile(A_birth_no, 33.3),  np.percentile(A_birth_no, 66.7),  np.amax(A_birth_no)]
    Abins = [[] for h in range(len(Asplit)-1)]
    dABins = [[] for h in range(len(Asplit)-1)]
    for h in range(len(A_birth_no)):
        for l in range(len(Asplit)-1):
            if (A_birth_no[h] >= Asplit[l] and A_birth_no[h] < Asplit[l + 1]): 
                Abins[l].append(A_birth_no[h]); dABins[l].append(A_div_no[h])   
    x = [np.median(Abins[h]) for h in range(len(Abins))]; y = [np.median(dABins[h]) for h in range(len(dABins))]
    yerr = [[y[h] - np.percentile(dABins[h], 25) for h in range(len(dABins))], [np.percentile(dABins[h], 75) - y[h] for h in range(len(dABins))]]     


    ax = fig.add_subplot(2, 3, 4)
    ax.scatter(A_birth_no, A_div_no, s = 100, color= 'c', marker ="o", edgecolors='None', alpha = 0.08) 
#ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
    m, b = np.polyfit(A_birth_no, A_div_no, 1); r, p = pearsonr(A_birth_no, A_div_no) 
    ax.errorbar(x, y, yerr =  yerr,  color = 'k', fmt='o')
    ax.plot(A_birth_no, m * np.array(A_birth_no) + b, linestyle = '-', color = 'r' )
   # lambda_, v = np.linalg.eig(np.cov(A_birth,A_div))
   # lambda_ = np.sqrt(lambda_) 
   # ell = Ellipse(xy=(np.mean(A_birth), np.mean(A_div)), width=lambda_[0]*2*excludeSTD_div, height=lambda_[1]*2*excludeSTD_div, angle=np.rad2deg(np.arccos(v[0, 0])))
   # ell.set_facecolor('none')
   # ax.add_artist(ell)
    leg = "n = " + str(len(A_birth_no)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "c.c., p = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
    #ax.text(0.07, 0.8, leg, fontsize=9, color = 'k', transform=ax.transAxes)
    #if(k ==3):
    #    ax.text(0.73, 0.84, '$(A_{\top})$', fontsize= 21, color = 'k', transform=ax.transAxes)
    #elif(k==5):
    #    ax.text(0.73, 0.84, '$(A_{\ta})$', fontsize= 21, color = 'k', transform=ax.transAxes)
    ax.set_xticks([-0.3, 0.0, 0.3])
    ax.set_yticks([-0.15, 0.0, 0.15])
    ax.set_xlim([-0.4, 0.4]); ax.set_ylim([-0.2, 0.2])
    ax.set_xlabel( '$\\alpha_b$', fontsize= 23)
    ax.set_ylabel( '$\\alpha_d$ ', fontsize= 23)
    #if(k ==3):
    #    ax.set_xlabel(  '$\\alpha_b  (A_{\top})$', fontsize= 23)
    #    ax.set_ylabel( '$\\alpha_d   (A_{\top})$ ', fontsize= 23)
    #elif(k ==5):
    #    ax.set_xlabel(  '$\\alpha_b  (A_{\ta})$', fontsize= 23)
    #    ax.set_ylabel( '$\\alpha_d   (A_{\ta})$ ', fontsize= 23)
    #else:
    #    ax.set_xlabel( '$\\alpha_b$', fontsize= 23)
    #    ax.set_ylabel( '$\\alpha_d$ ', fontsize= 23)


    A_birth_no = [A_birth[i] for i in range(len(A_birth)) if not outlier(A_birth, A_inc, A_birth[i], A_inc[i],excludeSTD_inc)]
    A_inc_no = [A_inc[i] for i in range(len(A_inc)) if not outlier(A_birth, A_inc, A_birth[i], A_inc[i], excludeSTD_inc)]
    all_plants_no = [all_plants[i] for i in range(len(all_plants)) if not outlier(A_birth, A_inc, A_birth[i], A_inc[i], excludeSTD_inc)]

    Abins = [[] for h in range(len(Asplit)-1)]
    dABins = [[] for h in range(len(Asplit)-1)]
    for h in range(len(A_birth_no)):
        for l in range(len(Asplit)-1):
            if (A_birth_no[h] >= Asplit[l] and A_birth_no[h] < Asplit[l + 1]): 
                Abins[l].append(A_birth_no[h]); dABins[l].append(A_inc_no[h])   
    x = [np.median(Abins[h]) for h in range(len(Abins))]; y = [np.median(dABins[h]) for h in range(len(dABins))]
    yerr = [[y[h] - np.percentile(dABins[h], 25) for h in range(len(dABins))], [np.percentile(dABins[h], 75) - y[h] for h in range(len(dABins))]]   

    ax = fig.add_subplot(2, 3, 5)
    ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
    ax.scatter(A_birth_no, A_inc_no, s = 100, color= 'c', marker ="o", edgecolors='None', alpha = 0.08) 
    #ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
    m, b = np.polyfit(A_birth_no, A_inc_no, 1); r, p = pearsonr(A_birth_no, A_inc_no) 
    ax.plot(A_birth_no, m * np.array(A_birth_no) + b, linestyle = '-', color = 'r' )

    leg = "n = " + str(len(A_birth_no)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "c.c., p = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
    #if(k ==3):
    #    ax.text(0.73, 0.84, '$(A_{\top})$', fontsize=21, color = 'k', transform=ax.transAxes)
    #elif(k ==5):
    #    ax.text(0.73, 0.84, '$(A_{\ta})$', fontsize=21, color = 'k', transform=ax.transAxes)
    #ax.text(0.07, 0.8, leg, fontsize=9, color = 'k', transform=ax.transAxes)
    ax.set_xticks([-0.3, 0.0, 0.3])
    ax.set_yticks([-0.3, 0.0, 0.3])
    ax.set_xlim([-0.4, 0.4]); ax.set_ylim([-0.4, 0.4])
    ax.set_xlabel(  '$\\alpha_b$', fontsize= 25)
    ax.set_ylabel( '$\\alpha_{\Delta}$ ', fontsize= 25)
    #if(k ==3):
    #    ax.set_xlabel(  '$\\alpha_b  (A_{\top})$', fontsize= 23)
    #    ax.set_ylabel( '$\\alpha_{\Delta}  (A_{\top})$ ', fontsize= 23)
    #elif(k ==5):
    #    ax.set_xlabel(  '$\\alpha_b  (A_{\ta})$', fontsize= 23)
    #    ax.set_ylabel( '$\\alpha_{\Delta}  (A_{\ta})$ ', fontsize= 23)
    #else:
    #    ax.set_xlabel( '$\\alpha_b$', fontsize= 23)
    #    ax.set_ylabel( '$\\alpha_{\Delta}$ ', fontsize= 23)

    A_birth_no = [A_birth[i] for i in range(len(A_birth)) if not outlier(A_birth, A_T, A_birth[i], A_T[i], excludeSTD_T)]
    A_T_no = [A_T[i] for i in range(len(A_T)) if not outlier(A_birth, A_T, A_birth[i], A_T[i], excludeSTD_T)]
    all_plants_no = [all_plants[i] for i in range(len(all_plants)) if not outlier(A_birth, A_T, A_birth[i], A_T[i], excludeSTD_T)]
    if k == 1: np.savetxt(resultsPath + "alphab_vs_alphaT_r_%s"%(radius) + ".out", [A_birth_no, A_T_no], delimiter=",", fmt="%s")


    Abins = [[] for h in range(len(Asplit)-1)]
    dABins = [[] for h in range(len(Asplit)-1)]
    for h in range(len(A_birth_no)):
        for l in range(len(Asplit)-1):
            if (A_birth_no[h] >= Asplit[l] and A_birth_no[h] < Asplit[l + 1]): 
                Abins[l].append(A_birth_no[h]); dABins[l].append(A_T_no[h])   
    x = [np.median(Abins[h]) for h in range(len(Abins))]; y = [np.median(dABins[h]) for h in range(len(dABins))]
    yerr = [[y[h] - np.percentile(dABins[h], 25) for h in range(len(dABins))], [np.percentile(dABins[h], 75) - y[h] for h in range(len(dABins))]]   
    ax = fig.add_subplot(2, 3, 6)
    ax.scatter(A_birth_no, A_T_no,  s = 100, color= 'c', marker ="o", edgecolors='None', alpha = 0.08) 
    #ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
    m, b = np.polyfit(A_birth_no, A_T_no, 1); r, p = pearsonr(A_birth_no, A_T_no) 
    ax.plot(A_birth_no, m * np.array(A_birth_no) + b, linestyle = '-', color = 'r' )
    ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
    leg = "n = " + str(len(A_birth_no)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "c.c., p = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
    #if(k ==3):
    #    ax.text(0.73, 0.84, '$(A_{\top})$', fontsize=21, color = 'k', transform=ax.transAxes)
    #elif(k ==5):
    #    ax.text(0.73, 0.84, '$(A_{\ta})$', fontsize=21, color = 'k', transform=ax.transAxes)
    #ax.text(0.07, 0.8, leg, fontsize=9, color = 'k', transform=ax.transAxes)
    #if(k ==3):
    #    ax.set_xlabel(  '$\\alpha_b  (A_{\top})$', fontsize= 23)
    #    ax.set_ylabel( '$\\alpha_T  (A_{\top})$ ', fontsize= 23)
    #elif(k ==5):
    #    ax.set_xlabel(  '$\\alpha_b  (A_{\ta})$', fontsize= 23)
    #    ax.set_ylabel( '$\\alpha_T  (A_{\ta})$ ', fontsize= 23)
    #else:
    #    ax.set_xlabel( '$\\alpha_b$', fontsize= 23)
    #    ax.set_ylabel( '$\\alpha_T$ ', fontsize= 23)
    ax.set_xlabel( '$\\alpha_b$', fontsize= 23)
    ax.set_ylabel( '$\\alpha_T$ ', fontsize= 23)
    ax.set_xticks([-0.3, 0.0, 0.3])
    ax.set_yticks([-0.3, 0.0, 0.3])
    ax.set_xlim([-0.4, 0.4]); ax.set_ylim([-0.4, 0.4])

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.48, hspace=None)
    fig.set_size_inches(13.0, 7.5)
    plt.savefig(resultsPath + "%sAsymBirth_%s_STD_%s.pdf"%(cell_vars[k], radius, excludeSTD_T), format='pdf', dpi=200)
    plt.close(fig)


    Sum_sis_birth_no = [Sum_sis_birth[i] for i in range(len(Sum_sis_birth)) if not outlier(Sum_sis_birth, Sum_sis_div, Sum_sis_birth[i], Sum_sis_div[i],excludeSTD_sumsis)]
    Sum_sis_div_no = [Sum_sis_div[i] for i in range(len(Sum_sis_div)) if not outlier(Sum_sis_birth, Sum_sis_div, Sum_sis_birth[i], Sum_sis_div[i],excludeSTD_sumsis)]
    all_plants_no = [all_plants[i] for i in range(len(all_plants)) if not outlier(Sum_sis_birth, Sum_sis_div, Sum_sis_birth[i], Sum_sis_div[i],excludeSTD_sumsis)]

    Asplit = [ np.amin(Sum_sis_birth_no), np.percentile(Sum_sis_birth_no, 33.3),  np.percentile(Sum_sis_birth_no, 66.7),  np.amax(Sum_sis_birth_no)]
    Abins = [[] for h in range(len(Asplit)-1)]
    dABins = [[] for h in range(len(Asplit)-1)]
    for h in range(len(Sum_sis_birth_no)):
        for l in range(len(Asplit)-1):
            if (Sum_sis_birth_no[h] >= Asplit[l] and Sum_sis_birth_no[h] < Asplit[l + 1]): 
                Abins[l].append(Sum_sis_birth_no[h]); dABins[l].append(Sum_sis_div_no[h])   
    x = [np.median(Abins[h]) for h in range(len(Abins))]; y = [np.median(dABins[h]) for h in range(len(dABins))]
    yerr = [[y[h] - np.percentile(dABins[h], 25) for h in range(len(dABins))], [np.percentile(dABins[h], 75) - y[h] for h in range(len(dABins))]]   

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(Sum_sis_birth_no, Sum_sis_div_no,  color="0.3",  marker ="o", edgecolors='None', alpha = 0.3) 
    #ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
    m, b = np.polyfit(Sum_sis_birth_no, Sum_sis_div_no, 1); r, p = pearsonr(Sum_sis_birth_no, Sum_sis_div_no) 
    ax.plot(Sum_sis_birth_no, m * np.array(Sum_sis_birth_no) + b, linestyle = '-', color = 'b' )
    ax.errorbar(x, y, yerr =  yerr, linestyle='-', color = 'k', marker='o')
    leg = "n = " + str(len(Sum_sis_birth_no)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "c.c., p = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
    ax.text(0.07, 0.8, leg, fontsize=8, color = 'k', transform=ax.transAxes)
    ax.set_xlabel(  '$(V_b^{\, sis1} + V_b^{\, sis2})/\mu_{V}$', fontsize= 12)
    ax.set_ylabel(  '$(V_d^{\, sis1} + V_d^{\, sis2})/\mu_{V}$', fontsize= 12)
    ax.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    ax.set_yticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
    ax.set_xlim([1., 2.5]); ax.set_ylim([1.5, 4.5])
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.28, hspace=None)

    fig.set_size_inches(6.0, 5.0)
    plt.savefig(resultsPath + "%s,SumSisBirth_Vs_sumSisDiv_%s_STD_%s.pdf"%(cell_vars[k], radius, excludeSTD_sumsis), format='pdf', dpi=200)
    plt.close(fig)


all_tdatd = all_time_of_day_at_division[0] + all_time_of_day_at_division[1] 
plt.hist(np.array(all_tdatd), bins =6)
plt.xlabel(' time of day at division (sunrise at 0hrs)')
plt.title(region)
plt.savefig(resultsPath + "times_of_day_at_division_%s.pdf"%radius, format='pdf', dpi=300)
plt.close()

all_tdatb = all_time_of_day_at_birth[0] + all_time_of_day_at_birth[1] 
plt.hist(np.array(all_tdatb), bins =6)
plt.xlabel(' time of day at division (sunrise at 0hrs)')
plt.title(region)
plt.savefig(resultsPath + "times_of_day_at_birth_%s.pdf"%radius, format='pdf', dpi=300)
plt.close()

plt.hist(np.array(all_tdatd), bins =3)
plt.xlabel(' time of day at division (sunrise at 0hrs)')
plt.title(region)
plt.savefig(resultsPath + "times_of_day_at_division_3bins_%s.pdf"%radius, format='pdf', dpi=300)
plt.close()


plt.hist(np.array(all_tdatb), bins =3)
plt.xlabel(' time of day at birth (sunrise at 0hrs)')
plt.title(region)
plt.savefig(resultsPath + "times_of_day_at_birth_3bins_%s.pdf"%radius, format='pdf', dpi=300)
plt.close()

pos = [4, 12, 20]
data, labelsPos = [], []
for p in pos:
    data.append([all_tdatd[i] - all_tdatb[i] for i in range(len(all_tdatd)) if all_tdatb[i] >= p - pos[0] and all_tdatb[i] < p + pos[0]])
    labelsPos.append(str(p - pos[0]) + "-" + str(p - pos[1]))

plt.boxplot(data, positions = pos, widths = 0.7)
plt.xlabel(' time of day at birth (sunrise at 0hrs)')
plt.ylabel(' time of day at division - time of day at birth')
plt.title(region)
plt.savefig(resultsPath + "times_day_birth_vs_time_day_div_min_birth_%s.pdf"%radius, format='pdf', dpi=300)
plt.close()
 
