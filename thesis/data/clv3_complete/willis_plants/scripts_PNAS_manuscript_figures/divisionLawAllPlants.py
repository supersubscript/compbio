import cPickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from scipy.stats.stats import skew
import all_plant_variables as apv 
from matplotlib.patches import Ellipse

# # # # # # # # # # # # # CELLS TO INCLUDE IN STATS: CENTRAL ZONE ONLY AND OUTLIERS
radius = '30'
region = '_initial central L1'     
excludeSTD = 2.5 #exclude cells outside this zone of standard deviation

# split data according to: "all_data", "early_times", "late_times", "inner_zone", "outer_zone"
resultsPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_individual_plant_data/all_data/"

# # # # # # # # # # # # # VARIABLE LABELS 
cell_vars = ["V", "V norm", "A", "Aop",  "Aip", "Aa"]

# # # # # # # # # # # # # FIGURES AND FIGURE VARIABLES
plantFontSize = 8
plt.locator_params(axis = 'x', nbins=4); plt.locator_params(axis = 'y', nbins=5)
plt.rcParams['xtick.labelsize'] = 7; plt.rcParams['ytick.labelsize'] = 7

figVarBVsVarD = [plt.figure(i) for i in range(len(cell_vars))]  
figVarBVsInc = [plt.figure(i + len(cell_vars)) for i in range(len(cell_vars))] 
figVarBVsVarDMod = [plt.figure(i + 2*len(cell_vars)) for i in range(len(cell_vars))]  
figVarBVsIncMod = [plt.figure(i + 3*len(cell_vars)) for i in range(len(cell_vars))]  

# # # # # # # # # # # # # FORMAT FOR EXPORTED DATA
varBVsVarD_fits = [["plant", " var bir. vs var div.", "slope", "intercept", "corr. coeff.", "p-value", " var bir. vs inc.", "slope", "intercept", "corr. coeff.", "p-value", "mean var birth", "coeff. var. var birth", "mean var div.", "std var division/mean var birth", "std (var S1 - var S2)/ 2(var S1 + var S2)"]]

varBVsVarD_fits_no_outliers = [["plant", " var bir. vs var div.", "slope", "intercept", "corr. coeff.", "p-value", " var bir. vs inc.", "slope", "intercept", "corr. coeff.", "p-value", "mean var birth", "coeff. var. var birth", "mean var div.", "std var division/mean var birth", "std (var S1 - var S2)/ 2(var S1 + var S2)"]]

def round_sig(x, sig=2):
    return np.round(x, sig-int(np.floor(np.log10(x)))-1)

def outlier(sampleInd, varInd):
    x0 = size_vars_birth[varInd][sampleInd]/mean_vars_birth[varInd] - np.mean(size_vars_birth[varInd]/mean_vars_birth[varInd]);
    y0 = size_vars_div[varInd][sampleInd]/mean_vars_birth[varInd] - np.mean(size_vars_div[varInd]/mean_vars_birth[varInd]);
    x1 =  np.cos(theta)*x0 - np.sin(theta)*y0; y1 =  np.sin(theta)*x0 + np.cos(theta)*y0
    if x1/lambda_[0] * x1/lambda_[0] + y1/lambda_[1] * y1/lambda_[1] < excludeSTD*excludeSTD: return False
    else : return True
    

for plant in range(len(apv.allPlants)):
    rootPath = apv.allRootPaths[plant]; timesList = apv.allTimesList[plant]; ALTDataPathList = apv.allALTDataPathList[plant]; segDataPathList = apv.allSegDataPathList[plant]
    meanVolumes = apv.allMeanVolumes[plant]; L1centre = apv.allL1centre[plant]

    dataFolder = rootPath
    times = [timesList[i][4] for i in range(0, len(timesList))]

    data_t0 = np.loadtxt(dataFolder + "Div_t0_new_" + radius + region + ".dat").astype(int) #(j0, m0, d0), cell divides between time[j0] with mother m0 and time[j0+1] with daughter d0
    data_t1 = np.loadtxt(dataFolder + "Div_t1_new_" + radius + region + ".dat").astype(int) #(j1, m1, d1), cell divides between time[j1] with mother m1 and time[j1+1] with daughter d1

    size_vars = [[cell_vars[i], []] for i in range(len(cell_vars))]
    size_vars_div, asymmetry_div = [[] for i in range(len(cell_vars))], [[] for i in range(len(cell_vars))]
    dist_tbirth, dist_centre = [], []
    time_since_birth, distance_over_cycle = [], []
    time_at_birth, time_of_day_at_birth  = [], []
    cell_cycle_time = []

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
                time_at_birth.append(times[j] - times[0])
                time_of_day_at_birth.append(times[j]%24)
                cell_cycle_time.append(times[data_t1[i][0] + 1] - times[j])
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
                #if(neighbourIndex < d0): sArea = sArea + dataSeg[0]['wall_surface'][(neighbourIndex, d0)]
                if(neighbourIndex > d0): sArea = sArea + dataSeg[0]['wall_surface'][(neighbourIndex, d0)]
                else: sArea = sArea + dataSeg[0]['wall_surface'][(d0, neighbourIndex)]
                if 1 not in dataSeg[0]['neigbourhood'][neighbourIndex] and neighbourIndex != 1: 
                    #if(neighbourIndex < d0): oPArea = oPArea + dataSeg[0]['wall_surface'][(neighbourIndex, d0)]
                    if(neighbourIndex > d0): oPArea = oPArea + dataSeg[0]['wall_surface'][(neighbourIndex, d0)]
                    else : oPArea = oPArea + dataSeg[0]['wall_surface'][(d0, neighbourIndex)]
            size_vars[2][1][-1].append(sArea); size_vars[3][1][-1].append(pArea); size_vars[4][1][-1].append(oPArea); size_vars[5][1][-1].append(sArea - pArea - oPArea)   

        #find tracked d0 at next timepoint and set equal to d0. Should be only one mother because of ways cells have been selected ... 
            fobj = file(ALTDataPathList[j])
            dataLineage, scores = cPickle.load(fobj)
            fobj.close()
            newDaughters = [x[0] for x in dataLineage if x[1] == d0]
            if (j != data_t1[i][0]) and (len(newDaughters) != 1) : print " More or less than one offspring ...  "
            if len(newDaughters) == 0 : print " Zero offspring ...  "; break
            else : d0 = newDaughters[0] 
        dist_centre.append(np.mean(dist_evolution))
        fobj = file(segDataPathList[data_t1[i][0] + 1])
        dataSeg = cPickle.load(fobj)
                
        vol, volNorm, sArea, pArea, oPArea = [0 for j in range(len(newDaughters))], [0 for j in range(len(newDaughters))], [0 for j in range(len(newDaughters))], [0 for j in range(len(newDaughters))], [0 for j in range(len(newDaughters))]

        for j in range(len(newDaughters)): 
            vol[j] = dataSeg[0]['volumes'][newDaughters[j]]; volNorm[j] = dataSeg[0]['volumes'][newDaughters[j]]/meanVolumes[data_t1[i][0] + 1]
            if 1 not in dataSeg[0]['neigbourhood'][newDaughters[j]] : print "No periclinal wall !!!!!!!!!!!!!!!!!!!!!! "
            #else: pArea[j] = dataSeg[0]['wall_surface'][(1, newDaughters[j])]
            else: pArea[j] = dataSeg[0]['wall_surface'][(newDaughters[j], 1)]
            for neighbourIndex in dataSeg[0]['neigbourhood'][newDaughters[j]]:
                #if(neighbourIndex < newDaughters[j]): sArea[j] = sArea[j] + dataSeg[0]['wall_surface'][(neighbourIndex, newDaughters[j])]
                if(neighbourIndex > newDaughters[j]): sArea[j] = sArea[j] + dataSeg[0]['wall_surface'][(neighbourIndex, newDaughters[j])]
                else: sArea[j] = sArea[j] + dataSeg[0]['wall_surface'][(newDaughters[j], neighbourIndex)]
                if 1 not in dataSeg[0]['neigbourhood'][neighbourIndex] and neighbourIndex != 1: 
                    #if(neighbourIndex < newDaughters[j]): oPArea[j] = oPArea[j] + dataSeg[0]['wall_surface'][(neighbourIndex, newDaughters[j])]
                    if(neighbourIndex > newDaughters[j]): oPArea[j] = oPArea[j] + dataSeg[0]['wall_surface'][(neighbourIndex, newDaughters[j])]
                    else : oPArea[j] = oPArea[j] + dataSeg[0]['wall_surface'][(newDaughters[j], neighbourIndex)]
        
        size_vars_div[0].append(np.sum(vol)); size_vars_div[1].append(np.sum(volNorm)); size_vars_div[2].append(np.sum(sArea)); size_vars_div[3].append(np.sum(pArea)); size_vars_div[4].append(np.sum(oPArea)); size_vars_div[5].append(np.sum(sArea) - np.sum(pArea) - np.sum(oPArea)) 

        if len(newDaughters) == 2 : 
            random = 1 if np.random.rand() < 0.5 else -1
            asymmetry_div[0].append(0.5*(vol[0] - vol[1])*random/np.sum(vol)); asymmetry_div[1].append(0.5*(volNorm[0] - volNorm[1])*random/np.sum(volNorm));  asymmetry_div[2].append(0.5*(sArea[0] - sArea[1])*random/np.sum(sArea));  asymmetry_div[4].append(0.5* random*(oPArea[0] -oPArea[1])/np.sum(oPArea)); asymmetry_div[5].append(0.5* random*(sArea[0] - sArea[1] - (pArea[0] - pArea[1]) -  (oPArea[0] - oPArea[1]))/(np.sum(sArea) - np.sum(pArea) - np.sum(oPArea))) 
            if np.sum(pArea) != 0 : asymmetry_div[3].append(0.5* random*(pArea[0] - pArea[1])/np.sum(pArea))
            else : asymmetry_div[3].append(0)
        else: 
            for asymmetry in asymmetry_div: asymmetry.append(0)
            print " Aberant number of daughters =  ", len(newDaughters)
   
    for var in size_vars: var[1] = np.array(var[1])
    size_vars_birth = []
    for var in size_vars: size_vars_birth.append(np.array([var[1][i][0] for i in range(0, len(var[1]))]))
    mean_vars_birth = []
    for var_birth in size_vars_birth: mean_vars_birth.append(np.mean(var_birth))
    mean_vars_birth[1] = 1.0
    mean_period = np.mean(cell_cycle_time)

    print " vs  ... "
    print "Var./mean at birth vs Var./mean at division ,  slope,  intercept,   corr. coeff., p-value  "
    for i in range(len(size_vars)):
        size_vars_birth[i]  =np.array(size_vars_birth[i]); size_vars_div[i] = np.array(size_vars_div[i])
        m, b = np.polyfit(size_vars_birth[i]/mean_vars_birth[i], size_vars_div[i]/mean_vars_birth[i], 1)
        r, p = pearsonr(size_vars_birth[i]/mean_vars_birth[i], size_vars_div[i]/mean_vars_birth[i])
        mI, bI = np.polyfit(size_vars_birth[i]/mean_vars_birth[i], np.array(size_vars_div[i]/mean_vars_birth[i]) - np.array(size_vars_birth[i]/mean_vars_birth[i]), 1)
        rI, pI = pearsonr(size_vars_birth[i]/mean_vars_birth[i],  np.array(size_vars_div[i]/mean_vars_birth[i]) - np.array(size_vars_birth[i]/mean_vars_birth[i]))
        print size_vars[i][0],",", m,",",  b,",", r, ",", p 
        varBVsVarD_fits.append([apv.allPlants[plant], size_vars[i][0], m, b, r, p, size_vars[i][0], mI, bI, rI, pI, np.mean(size_vars_birth[i]), np.std(size_vars_birth[i]/mean_vars_birth[i]), np.mean(size_vars_div[i]), np.std(size_vars_div[i]/mean_vars_birth[i]), np.std(asymmetry_div[i])])
        ax = figVarBVsVarD[i].add_subplot(2,3, plant + 1)
        ax.scatter(size_vars_birth[i]/mean_vars_birth[i], size_vars_div[i]/mean_vars_birth[i], c = dist_centre, marker ="o", edgecolors='None', alpha = 0.3)
        ax.plot(size_vars_birth[i]/mean_vars_birth[i], m*size_vars_birth[i]/mean_vars_birth[i] + b, '--k')
        x = size_vars_birth[i]/mean_vars_birth[i]; y = size_vars_div[i]/mean_vars_birth[i]
        lambda_, v = np.linalg.eig(np.cov(x,y))
        lambda_ = np.sqrt(lambda_) 
        ell = Ellipse(xy=(np.mean(x), np.mean(y)), width=lambda_[0]*2*excludeSTD, height=lambda_[1]*2*excludeSTD, angle=np.rad2deg(np.arccos(v[0, 0])))
        ell.set_facecolor('none')
        ax.add_artist(ell)
       # ax.set_xlim(axLimitsVarVsCorrT[i][0]); ax.set_ylim(axLimitsVarVsCorrT[i][1])
        leg = "n = " + str(len(size_vars_birth[i])) + "\n"  + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "c.c., p = "  + str(round(r, 2)) + ", " + str(round_sig(p)) + "\n"
        ax.text(0.07, 0.79, leg, fontsize=7, color = 'k', transform=ax.transAxes)
        ax.set_title(apv.allPlants[plant], fontsize = plantFontSize)

        ax = figVarBVsInc[i].add_subplot(2,3, plant + 1)
        ax.scatter(size_vars_birth[i]/mean_vars_birth[i], size_vars_div[i]/mean_vars_birth[i] - size_vars_birth[i]/mean_vars_birth[i], c = dist_centre, marker ="o", edgecolors='None', alpha = 0.3)
        ax.plot(size_vars_birth[i]/mean_vars_birth[i], mI*size_vars_birth[i]/mean_vars_birth[i] + bI, '--k')
       # ax.set_xlim(axLimitsVarVsCorrT[i][0]); ax.set_ylim(axLimitsVarVsCorrT[i][1])
        leg = "n = " + str(len(size_vars_birth[i])) + "\n"  + "sl., int. = "  + str(round(mI,2)) + ", " + str(round(bI,2)) + "\n" + "c.c., p = "  + str(round(rI, 2)) + ", " + str(round_sig(pI)) + "\n"
        ax.text(0.07, 0.79, leg, fontsize=7, color = 'k', transform=ax.transAxes)
        ax.set_title(apv.allPlants[plant], fontsize = plantFontSize)

    for i in range(len(size_vars)):
        x = size_vars_birth[i]/mean_vars_birth[i]; y = size_vars_div[i]/mean_vars_birth[i]
        lambda_, v = np.linalg.eig(np.cov(x,y))
        lambda_ = np.sqrt(lambda_) 
        ell = Ellipse(xy=(np.mean(x), np.mean(y)), width=lambda_[0]*2*excludeSTD, height=lambda_[1]*2*excludeSTD, angle=np.rad2deg(np.arccos(v[0, 0])))
        ell.set_facecolor('none')

        theta = -np.arccos(v[0, 0])
        sBMod = []; sDMod = []; asymMod = []; distMod = []; timeBirthMod = []
        for j in range(len(size_vars_birth[i])):
            if not outlier(j, i):
                sBMod.append(size_vars_birth[i][j]); sDMod.append(size_vars_div[i][j]); asymMod.append(asymmetry_div[i][j]); distMod.append(dist_centre[j]); timeBirthMod.append(time_at_birth[j])
        sBMod = np.array(sBMod); sDMod = np.array(sDMod)
        print " fraction included : ", len(sBMod)/float(len(size_vars_birth[i]))  
        mSB = np.mean(sBMod) if i != 1 else 1.0
        m, b = np.polyfit(sBMod/mSB, sDMod/mSB, 1)
        r, p = pearsonr(sBMod/mSB, sDMod/np.mean(sBMod))
        mI, bI = np.polyfit(sBMod/mSB, sDMod/mSB - sBMod/mSB, 1)
        rI, pI = pearsonr(sBMod/mSB, sDMod/mSB - sBMod/mSB)
        print size_vars[i][0],",", m,",",  b,",", r, ",", p 
        varBVsVarD_fits_no_outliers.append([apv.allPlants[plant], size_vars[i][0], m, b, r, p, size_vars[i][0], mI, bI, rI, pI, np.mean(sBMod), np.std(sBMod/mSB), np.mean(sDMod), np.std(sDMod/mSB), np.std(asymMod)])
        ax = figVarBVsVarDMod[i].add_subplot(2,3, plant + 1)
        ax.scatter(sBMod/mean_vars_birth[i], sDMod/mean_vars_birth[i], c = distMod, marker ="o", edgecolors='None', alpha = 0.3)
        ax.plot(sBMod/mean_vars_birth[i], m*sBMod/mean_vars_birth[i] + b, '--k')
        ax.add_artist(ell)
       # ax.set_xlim(axLimitsVarVsCorrT[i][0]); ax.set_ylim(axLimitsVarVsCorrT[i][1])
        leg = "n = " + str(len(sBMod)) + "\n"  + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "c.c., p = "  + str(round(r, 2)) + ", " + str(round_sig(p)) + "\n"
        ax.text(0.07, 0.75, leg, fontsize=7, color = 'k', transform=ax.transAxes)
        ax.set_title(apv.allPlants[plant], fontsize = plantFontSize)

        ax = figVarBVsIncMod[i].add_subplot(2,3, plant + 1)
        ax.scatter(sBMod/mSB, sDMod/mSB - sBMod/mSB, c = distMod, marker ="o", edgecolors='None', alpha = 0.3)
        ax.plot(np.array(sBMod)/mSB, mI*np.array(sBMod)/mSB + bI, '--k')
       # ax.set_xlim(axLimitsVarVsCorrT[i][0]); ax.set_ylim(axLimitsVarVsCorrT[i][1])
        leg = "n = " + str(len(sBMod)) + "\n"  + "sl., int. = "  + str(round(mI,2)) + ", " + str(round(bI,2)) + "\n" + "c.c., p = "  + str(round(rI, 2)) + ", " + str(round_sig(pI)) + "\n"
        ax.text(0.07, 0.75, leg, fontsize=7, color = 'k', transform=ax.transAxes)
        ax.set_title(apv.allPlants[plant], fontsize = plantFontSize)


#sort data to lump cell parameters together
varBVsVarD_fits_sort = [varBVsVarD_fits[0]]
for entry in size_vars:
    for row in varBVsVarD_fits : 
        if entry[0] == row[1] : varBVsVarD_fits_sort.append(row) 
np.savetxt(resultsPath + "varBVsVarD_fits.out", np.array(varBVsVarD_fits_sort), delimiter="\t", fmt="%s")

varBVsVarD_fits_no_outliers_sort = [varBVsVarD_fits_no_outliers[0]]
for entry in size_vars:
    for row in varBVsVarD_fits_no_outliers : 
        if entry[0] == row[1] : varBVsVarD_fits_no_outliers_sort.append(row) 
np.savetxt(resultsPath + "varBVsVarD_no_oultiers_fits.out", np.array(varBVsVarD_fits_no_outliers_sort), delimiter="\t", fmt="%s")

for i in range(len(cell_vars)):
    plt.figure(i)
    figVarBVsVarD[i].text(0.5, 0.04, '%s at birth '%(cell_vars[i]), fontsize = 8,  ha='center', va='center')
    figVarBVsVarD[i].text(0.06, 0.5, '%s at division '%(cell_vars[i]),fontsize = 8,  ha='center', va='center', rotation='vertical')
    plt.savefig(resultsPath + "%sBirthVs%sDiv_%s.png"%(cell_vars[i], cell_vars[i], region), format='png', dpi=200)
    plt.close(i)
    plt.figure(i+len(cell_vars))
    figVarBVsInc[i].text(0.5, 0.04, '%s at birth '%(cell_vars[i]), fontsize = 8,  ha='center', va='center')
    figVarBVsInc[i].text(0.06, 0.5, '%s increment '%(cell_vars[i]),fontsize = 8,  ha='center', va='center', rotation='vertical')
    plt.savefig(resultsPath + "%sBirthVs%sInc_%s.png"%(cell_vars[i], cell_vars[i], region), format='png', dpi=200)
    plt.close(i + len(cell_vars))
    plt.figure(i + 2 * len(cell_vars))
    figVarBVsVarDMod[i].text(0.5, 0.04, '%s at birth '%(cell_vars[i]), fontsize = 8,  ha='center', va='center')
    figVarBVsVarDMod[i].text(0.06, 0.5, '%s at division '%(cell_vars[i]),fontsize = 8,  ha='center', va='center', rotation='vertical')
    plt.savefig(resultsPath + "%sBirthVs%sDiv_no_outliers_%s.png"%(cell_vars[i], cell_vars[i], region), format='png', dpi=200)
    plt.close(i + 2 * len(cell_vars))
    plt.close(i + len(cell_vars))
    plt.figure(i + 3 * len(cell_vars))
    figVarBVsIncMod[i].text(0.5, 0.04, '%s at birth  '%(cell_vars[i]), fontsize = 8,  ha='center', va='center')
    figVarBVsIncMod[i].text(0.06, 0.5, '%s increment '%(cell_vars[i]),fontsize = 8,  ha='center', va='center', rotation='vertical')
    plt.savefig(resultsPath + "%sBirthVs%sInc_no_outliers_%s.png"%(cell_vars[i], cell_vars[i], region), format='png', dpi=200)
    plt.close(i+ 3 * len(cell_vars))
