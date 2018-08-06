import cPickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.stats import skew, skewtest
import statsmodels.api as sm

##################################################
from config_2015_02_16_plant1 import ALTDataPathList, segDataPathList, timesList, savePath
from indices_parameters_2015_02_16_plant1 import L1centre
radius = 30
zMin = 8
region = 'L1'
label = 'p1'
t_begin = 0
spTOfDay = [0,4, 8,12, 16,20, 24] 
################################################
plt.rcParams['xtick.labelsize'] = 14; plt.rcParams['ytick.labelsize'] = 14


times = [timesList[i][4] for i in range(0, len(timesList))]
print times
timesLabel = [float(timesList[i][0][0:-3]) for i in range(0, len(timesList))]
timesLabel2 = [timesList[i][0][0:-3] for i in range(0, len(timesList))]

for i in range(0, int(len(timesLabel2)/2)): timesLabel2[2*i + 1] = ''
print timesLabel
print timesLabel2

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
            if(neighbourIndex > cellInd): ipArea = ipArea + segData[tInd][0]['wall_surface'][(neighbourIndex, cellInd)]
            else : ipArea = ipArea + segData[tInd][0]['wall_surface'][(cellInd, neighbourIndex)]
    return ipArea


def wallArea(tInd, cellInd):
    sArea = 0
    for neighbourIndex in segData[tInd][0]['neigbourhood'][cellInd]:
        if(neighbourIndex > cellInd): sArea = sArea + segData[tInd][0]['wall_surface'][(neighbourIndex, cellInd)]
        else: sArea = sArea + segData[tInd][0]['wall_surface'][(cellInd, neighbourIndex)]
    return sArea

def round_sig(x, sig=2):
    if x ==0 : return 0.0
    elif x > 0: return np.round(x, sig-int(np.floor(np.log10(x)))-1)
    else: return -1*np.round(-x, sig-int(np.floor(np.log10(-x)))-1)


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


all_vols = []
all_num_L1neighs = []
all_opA = []
all_ipA = []
all_A = []
all_aA = []
all_dist = []
all_times = []
vols_t = []
numberDivisions = []
numberDivisionsPerCell = []
numberDivTofDay = [[] for x in range(len(spTOfDay)-1)]
volTofDay = [[] for x in range(len(spTOfDay)-1)]

for i in range(0, len(segDataPathList)):
    vols_t.append([])
    if i < len(segDataPathList) - 1: 
        #identify all L1 cells that are tracked to the next timepoint
        indices_dup = [x[1] for x in dataLineages[i] if 1 in segData[i][0]['neigbourhood'][x[1]]]
        indices = [indices_dup[j] for j in range(len(indices_dup)) if indices_dup[j] not in indices_dup[j+1:]]
        numberDivisions.append(0)   
    else: indices = [x[0] for x in dataLineages[i-1] if 1 in segData[i][0]['neigbourhood'][x[0]]]

    cellsInRange = []
    for index in indices: 
        xSq = (segData[i][0]['barycenter'][index][0] - L1centre[i][0]) * (segData[i][0]['barycenter'][index][0] - L1centre[i][0])
        ySq = (segData[i][0]['barycenter'][index][1] - L1centre[i][1]) * (segData[i][0]['barycenter'][index][1] - L1centre[i][1])
        zSq = (segData[i][0]['barycenter'][index][2] - L1centre[i][2]) * (segData[i][0]['barycenter'][index][2] - L1centre[i][2])
        if( np.sqrt(xSq + ySq + zSq) < radius and segData[i][0]['barycenter'][index][2] > zMin and times[i] > t_begin):             
            d = index         
            if(innPeriWallArea(i, d) != 0 and segData[i][0]['wall_surface'][(d,1)] != 0):
                all_num_L1neighs.append(len([x for x in segData[i][0]['neigbourhood'][d] if 1 in segData[i][0]['neigbourhood'][x]]))
                all_dist.append(np.sqrt(xSq + ySq + zSq))
                all_vols.append(segData[i ][0]['volumes'][d])
                all_opA.append(segData[i][0]['wall_surface'][(d,1)])
                all_A.append(wallArea(i, d))
                all_ipA.append(innPeriWallArea(i, d))
                all_aA.append(wallArea(i, d) -innPeriWallArea(i, d) -segData[i][0]['wall_surface'][(d,1)])
                all_times.append(times[i] - times[0])
                vols_t[-1].append(segData[i][0]['volumes'][d])
            if i < len(segDataPathList) - 1: 
                cellsInRange.append(index)
                numDaughters = [x[0] for x in dataLineages[i] if x[1] == index]
                if(len(numDaughters) > 1): 
                    numberDivisions[i] =  numberDivisions[i] +  len(numDaughters) 
    if i < len(segDataPathList) - 1: 
        numberDivisionsPerCell.append(float(numberDivisions[i])/len(cellsInRange))  
                
    print len(cellsInRange)
    for j in range(len(spTOfDay)-1):  
        if(i < len(segDataPathList) - 1 and times[i]%24 >= spTOfDay[j] and times[i]%24 < spTOfDay[j+1]):
            numberDivTofDay[j].append(numberDivisions[i])
            for vol in vols_t[-1]: volTofDay[j].append(vol)


sNumDivTofDay = [np.sum(numbDiv) for numbDiv in numberDivTofDay]

fig =plt.figure()   
ax = fig.add_subplot(1,1,1)
leg = "$N$ = " + str(np.sum(sNumDivTofDay)) + "\n" 
plt.bar(np.array(spTOfDay)[:-1], sNumDivTofDay, 0.35, alpha=0.5, color = 'c')
#plt.xlim(0, 2.2)
#plt.xlabel('$V/\mu_V$')
#plt.ylabel('frequency')
#adjustFigAspect(fig,aspect=0.8)
ax.text(0.07, 0.85, leg, fontsize=6, color = 'k', transform=ax.transAxes)
fig.set_size_inches(5.5, 4.0)
plt.gcf().subplots_adjust(bottom=0.12, left = 0.12)
#plt.yticks(np.arange(0.0, 401, 100))
plt.savefig(savePath + "hist_timeOfDayVsNumDivs_%s_z_%s_tb_%s_%s_fixed_range.pdf"%(radius, zMin, t_begin, region), format='pdf', dpi=300)
plt.close(fig)

labs = []
for t in range(len(spTOfDay)-1): labs.append(str(spTOfDay[t]) + '~' + str(spTOfDay[t+1]))

totS = np.sum([len(x) for x in volTofDay])
leg = "$N$ = " + str(totS) + "\n" 
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
sm.graphics.violinplot(volTofDay, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs', 'violin_lw':0,
                                  'label_fontsize':'large',
                                  'label_rotation':30})
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = (0, 4.5, 0, 800))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest', alpha =0.1, extent = (4.5, 7.6, 0, 800))
ax.set_xlabel(' phase of 24 h light/dark cycle (h)',  fontsize = 16)
#ax.set_ylabel(' $d%s/dt * T/\Delta_{%s}$ '%(cell_vars[h], cell_vars[h]))
ax.set_ylabel(' vol. ($\mu$m$^3$)', fontsize = 18)
ax.set_xlim([0, len(spTOfDay)])
ax.set_ylim([25,310])
ax.text(0.72, 0.82, leg, fontsize=15, color = 'k', transform=ax.transAxes)
ax.set_aspect('auto')
plt.gcf().subplots_adjust(bottom=0.24, left = 0.14)
fig.set_size_inches(6.5, 4.2)
plt.savefig(savePath + "%s_violin_timeOfDayVsVols_%s_z_%s_tb_%s_%s_fixed_range.pdf"%(label, radius, zMin, t_begin, region), format='pdf', dpi=300)
plt.close(fig)


colours = (np.array(all_dist) - np.min(all_dist))/(np.max(all_dist) - np.min(all_dist))
meanV =  np.mean(all_vols)
meanopA = np.mean(all_opA)
meanipA = np.mean(all_ipA)
meanA = np.mean(all_A)
meanaA = np.mean(all_aA)

print "mean V: ", meanV, " COV V : ", np.std(all_vols)/meanV, "skewness, p: ", skew(np.array(all_vols)/meanV), skewtest(np.array(all_vols)/meanV)

print "mean Aop: ", meanopA


Asplit = [ np.amin(np.log(np.array(all_vols)/meanV)), -0.2, 0.0, 0.2, np.amax(np.log(np.array(all_vols)/meanV))]
Abins = [[] for k in range(len(Asplit)-1)]; dABins = [[] for k in range(len(Asplit)-1)]
for k in range(len(all_vols)):
    for l in range(len(Asplit)-1):
        if (np.log(all_vols[k]/meanV) >= Asplit[l] and np.log(all_vols[k]/meanV) < Asplit[l + 1]): 
            Abins[l].append(np.log(all_vols[k]/meanV)); dABins[l].append(np.log(all_opA[k]/meanopA))   
x = [np.median(Abins[k]) for k in range(len(Abins))]; y = [np.median(dABins[k]) for k in range(len(dABins))]
yerr = [[y[k] - np.percentile(dABins[k], 25) for k in range(len(dABins))], [np.percentile(dABins[k], 75) - y[k] for k in range(len(dABins))]]           


m, b, r, p, std_err = linregress(np.log(np.array(all_vols)/meanV), np.log(np.array(all_opA)/meanopA))
leg = "$N$ = " + str(len(all_vols)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, p = " + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"

fig =plt.figure()   
ax = fig.add_subplot(1,1,1)
ax.scatter(np.log(np.array(all_vols)/meanV), np.log(np.array(all_opA)/meanopA),  c = 'c', marker ="o", edgecolors='None', alpha = 0.08, s= 40)
#plt.(0.0, 4.0)
#ax.plot(np.log(np.array(all_vols)/meanV), 0.9 * np.log(np.array(all_vols)/meanV) + 0, 'k-.')
ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt ='o')
ax.plot(np.arange(-.6, 0.9, 0.6), m * np.arange(-.6, 0.9, 0.6) + b, 'r-')
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(4.8, 4.2)
ax.text(0.07, 0.8, leg, fontsize=8, color = 'k', transform=ax.transAxes)
adjustFigAspect(fig,aspect=1.0)
plt.gcf().subplots_adjust(bottom=0.18, left = 0.28)
plt.ylabel('$\log( A_{\top}/\mu_{A_{\top}})$', fontsize = 19)
plt.xlabel('$\log( V/\mu_V)$', fontsize = 19)
plt.ylim(-0.9,0.8)
plt.xlim(-0.9,0.8)
plt.yticks(np.arange(-.6, 0.8, 0.4))
plt.xticks(np.arange(-.6, 0.8, 0.4))
plt.savefig(savePath + "logOPA_vs_logV_r_%s_z_%s_%s_fixed_range.pdf"%(radius, zMin, region), format='pdf', dpi=300)
plt.close(fig)


dABins = [[] for k in range(len(Asplit)-1)]
for k in range(len(all_vols)):
    for l in range(len(Asplit)-1):
        if (np.log(all_vols[k]/meanV) >= Asplit[l] and np.log(all_vols[k]/meanV) < Asplit[l + 1]): 
            dABins[l].append(np.log(all_A[k]/meanA))   
x = [np.median(Abins[k]) for k in range(len(Abins))]; y = [np.median(dABins[k]) for k in range(len(dABins))]
yerr = [[y[k] - np.percentile(dABins[k], 25) for k in range(len(dABins))], [np.percentile(dABins[k], 75) - y[k] for k in range(len(dABins))]]   

m, b, r, p, std_err = linregress(np.log(np.array(all_vols)/meanV), np.log(np.array(all_A)/meanA))
leg = "$N$ = " + str(len(all_vols)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, p = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"

fig =plt.figure()   
ax = fig.add_subplot(1,1,1)
ax.scatter(np.log(np.array(all_vols)/meanV), np.log(np.array(all_A)/meanA),  c = 'c', alpha = 0.08, s= 40, marker ="o", edgecolors='None')
#plt.(0.0, 4.0)
#ax.plot(np.log(np.array(all_vols)/meanV),m * np.log(np.array(all_vols)/meanV) + b, 'b-.')
ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt ='o')
ax.plot(np.arange(-.6, 0.9, 0.6), m * np.arange(-.6, 0.9, 0.6) + b, 'r-')
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(4.8, 4.2)
ax.text(0.07, 0.8, leg, fontsize=8, color = 'k', transform=ax.transAxes)
adjustFigAspect(fig,aspect=1.0)
plt.gcf().subplots_adjust(bottom=0.18, left = 0.28)
plt.ylabel('$\log( A/\mu_A)$', fontsize =19)
plt.xlabel('$\log( V/\mu_V)$', fontsize = 19)
plt.ylim(-0.9,0.8)
plt.xlim(-0.9,0.8)
plt.yticks(np.arange(-.6, 0.8, 0.4))
plt.xticks(np.arange(-.6, 0.8, 0.4))
plt.savefig(savePath + "logA_vs_logV_r_%s_z_%s_%s_fixed_range.pdf"%(radius, zMin, region), format='pdf', dpi=300)
plt.close()

dABins = [[] for k in range(len(Asplit)-1)]
for k in range(len(all_vols)):
    for l in range(len(Asplit)-1):
        if (np.log(all_vols[k]/meanV) >= Asplit[l] and np.log(all_vols[k]/meanV) < Asplit[l + 1]): 
            dABins[l].append(np.log(all_aA[k]/meanaA))   
x = [np.median(Abins[k]) for k in range(len(Abins))]; y = [np.median(dABins[k]) for k in range(len(dABins))]
yerr = [[y[k] - np.percentile(dABins[k], 25) for k in range(len(dABins))], [np.percentile(dABins[k], 75) - y[k] for k in range(len(dABins))]]   


m, b, r, p, std_err = linregress(np.log(np.array(all_vols)/meanV), np.log(np.array(all_aA)/meanaA))
leg = "$N$ = " + str(len(all_vols)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, p = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"

fig =plt.figure()   
ax = fig.add_subplot(1,1,1)
ax.scatter(np.log(np.array(all_vols)/meanV), np.log(np.array(all_aA)/meanaA),  c = 'c', alpha = 0.08, s =40, marker ="o", edgecolors='None')
#plt.(0.0, 4.0)
ax.plot(np.arange(-.6, 0.9, 0.6), m * np.arange(-.6, 0.9, 0.6) + b, 'r-')
ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(4.8, 4.2)
ax.text(0.07, 0.8, leg, fontsize=8, color = 'k', transform=ax.transAxes)
adjustFigAspect(fig,aspect=1.0)
plt.gcf().subplots_adjust(bottom=0.18, left = 0.28)
plt.ylabel('$\log( A_{\ta}/\mu_{A_{\ta}})$', fontsize = 19)
plt.xlabel('$\log( V/\mu_V)$', fontsize = 19)
plt.ylim(-0.45,0.45)
plt.xlim(-0.9,0.9)
plt.yticks(np.arange(-.3, 0.5, 0.3))
plt.xticks(np.arange(-.6, 0.9, 0.6))
plt.savefig(savePath + "logaA_vs_logV_r_%s_z_%s_%s_fixed_range.pdf"%(radius, zMin, region), format='pdf', dpi=300)
plt.close()




dABins = [[] for k in range(len(Asplit)-1)]
for k in range(len(all_vols)):
    for l in range(len(Asplit)-1):
        if (np.log(all_vols[k]/meanV) >= Asplit[l] and np.log(all_vols[k]/meanV) < Asplit[l + 1]): 
            dABins[l].append(np.log(all_ipA[k]/meanipA))   
x = [np.median(Abins[k]) for k in range(len(Abins))]; y = [np.median(dABins[k]) for k in range(len(dABins))]
yerr = [[y[k] - np.percentile(dABins[k], 25) for k in range(len(dABins))], [np.percentile(dABins[k], 75) - y[k] for k in range(len(dABins))]]   


m, b, r, p, std_err = linregress(np.log(np.array(all_vols)/meanV),np.log(np.array(all_ipA)/meanipA))
print m, b, r, p
leg = "$N$ = " + str(len(all_vols)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, p = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"

fig =plt.figure()   
ax = fig.add_subplot(1,1,1)
ax.scatter(np.log(np.array(all_vols)/meanV), np.log(np.array(all_ipA)/meanipA),  c= 'c', alpha = 0.08,  s = 40, marker ="o", edgecolors='None')
#plt.(0.0, 4.0)
#ax.plot(np.log(np.array(all_vols)/meanV), m * np.log(np.array(all_vols)/meanV) + b, 'b-.')
ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt ='o')
ax.plot(np.arange(-.6, 0.9, 0.6), m * np.arange(-.6, 0.9, 0.6) + b, 'r-')
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(4.8, 4.2)
ax.text(0.07, 0.8, leg, fontsize=8, color = 'k', transform=ax.transAxes)
adjustFigAspect(fig,aspect=1.0)
plt.gcf().subplots_adjust(bottom=0.18, left = 0.28)
plt.ylabel('$\log( A_{\tip}/\mu_{A_{\tip}})$', fontsize = 19)
plt.xlabel('$\log( V/\mu_V)$', fontsize = 19)
plt.ylim(-0.9,0.8)
plt.xlim(-0.9,0.8)
plt.yticks(np.arange(-.6, 0.8, 0.4))
plt.xticks(np.arange(-.6, 0.8, 0.4))
plt.savefig(savePath + "logipA_vs_logV_r_%s_z_%s_%s_fixed_range.pdf"%(radius, zMin, region), format='pdf', dpi=300)
plt.close()


colours = (np.array(all_vols) - np.min(all_vols))/(np.max(all_vols) - np.min(all_vols))
m, b, r, p, std_err = linregress(all_dist, np.array(all_vols)/meanV)
print m, b, r, p
leg = "$N$ = " + str(len(all_vols)) + "\n" + "sl., int. = "  + str(round_sig(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
fig =plt.figure()   
ax = fig.add_subplot(1,1,1)
ax.scatter(all_dist, np.array(all_vols)/meanV, c = 'c', alpha = 0.08,  s = 40, marker ="o", edgecolors='None')
#plt.(0.0, 4.0)
#ax.plot(all_dist, m * np.array(all_dist) + b, 'k-.')
plt.xlabel('distance from $O$ ($\mu$m)', fontsize = 18)
plt.ylabel('vol./mean(vol.)', fontsize = 18)
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(7.0, 5.3)
#ax.text(0.07, 0.77, leg, fontsize=13, color = 'k', transform=ax.transAxes)
adjustFigAspect(fig,aspect=1.0)
plt.gcf().subplots_adjust(bottom=0.19, left = 0.15)
plt.ylim(0.1,2.2)
plt.xlim(0, radius + 2)
plt.yticks(np.arange(0.5, 2.5, 0.5))
plt.xticks(np.arange(10, radius + 1, 10))
plt.savefig(savePath + "%s_dist_vs_V_r_%s_z_%s_%s_fixed_range.pdf"%(label, radius, zMin, region), format='pdf', dpi=300)
plt.close()

m, b, r, p, std_err = linregress(all_dist, np.array(all_A)/meanA)
print m, b, r, p
leg = "$N$ = " + str(len(all_A)) + "\n" + "sl., int. = "  + str(round_sig(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
fig =plt.figure()   
ax = fig.add_subplot(1,1,1)
ax.scatter(all_dist, np.array(all_A)/meanA, c= colours, alpha = 0.3, s=50, marker ="o", edgecolors='None')
plt.xlabel('distance from $O$ ($\mu$m)', fontsize = 18)
plt.ylabel('tot. area/mean(tot. area)', fontsize = 18)
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(7.0, 5.3)
ax.text(0.07, 0.77, leg, fontsize=13, color = 'k', transform=ax.transAxes)
adjustFigAspect(fig,aspect=1.0)
plt.gcf().subplots_adjust(bottom=0.19, left = 0.15)
plt.ylim(0.1,2.6)
plt.xlim(0, radius + 2)
plt.yticks(np.arange(0.5, 3.0, 0.5))
plt.xticks(np.arange(10, radius + 1, 10))
plt.savefig(savePath + "%s_dist_vs_A_r_%s_z_%s_%s_fixed_range.pdf"%(label, radius, zMin, region), format='pdf', dpi=300)
plt.close()


m, b, r, p, std_err = linregress(all_dist, np.array(all_opA)/meanopA)
print m, b, r, p
leg = "$N$ = " + str(len(all_vols)) + "\n" + "sl., int. = "  + str(round_sig(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"

fig =plt.figure()   
ax = fig.add_subplot(1,1,1)
ax.scatter(all_dist, np.array(all_opA)/meanopA, c= colours, alpha = 0.3,  marker ="o", s =50, edgecolors='None')
plt.xlabel('distance from $O$ ($\mu$m)', fontsize = 18)
plt.ylabel('outer peri. area/mean(outer peri. area)', fontsize = 18)
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(7.0, 5.3)
#ax.text(0.07, 0.77, leg, fontsize=13, color = 'k', transform=ax.transAxes)
adjustFigAspect(fig,aspect=1.0)
plt.gcf().subplots_adjust(bottom=0.19, left = 0.15)
plt.ylim(0.1,2.2)
plt.xlim(0, radius + 2)
plt.yticks(np.arange(0.5, 2.5, 0.5))
plt.xticks(np.arange(10, radius + 1, 10))
plt.savefig(savePath + "%s_dist_vs_opA_r_%s_z_%s_%s_fixed_range.pdf"%(label, radius, zMin, region), format='pdf', dpi=300)
plt.close()

m, b, r, p, std_err = linregress(all_dist, np.array(all_ipA)/meanipA)
print m, b, r, p
leg = "$N$ = " + str(len(all_vols)) + "\n" + "sl., int. = "  + str(round_sig(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"

fig =plt.figure()   
ax = fig.add_subplot(1,1,1)
ax.scatter(all_dist, np.array(all_ipA)/meanipA, c= colours, alpha = 0.3, s=50, marker ="o", edgecolors='None')
plt.xlabel('distance from $O$ ($\mu$m)', fontsize = 18)
plt.ylabel('inner peri. area/mean(outer peri. area)', fontsize = 18)
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(7.0, 5.3)
ax.text(0.07, 0.77, leg, fontsize=13, color = 'k', transform=ax.transAxes)
adjustFigAspect(fig,aspect=1.0)
plt.gcf().subplots_adjust(bottom=0.19, left = 0.15)
plt.ylim(0.1,2.6)
plt.xlim(0, radius + 2)
plt.yticks(np.arange(0.5, 3.0, 0.5))
plt.xticks(np.arange(10, radius + 1, 10))
plt.savefig(savePath + "%s_dist_vs_ipA_r_%s_z_%s_%s_fixed_range.pdf"%(label, radius, zMin, region), format='pdf', dpi=300)
plt.close()


m, b, r, p, std_err = linregress(all_dist, np.array(all_aA)/meanaA)
print m, b, r, p
leg = "$N$ = " + str(len(all_vols)) + "\n" + "sl., int. = "  + str(round_sig(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
fig =plt.figure()   
ax = fig.add_subplot(1,1,1)
ax.scatter(all_dist, np.array(all_aA)/meanaA, c= colours, alpha = 0.3, s =50, marker ="o", edgecolors='None')
plt.xlabel('distance from $O$ ($\mu$m)', fontsize = 18)
plt.ylabel('anti. peri. area/mean(outer peri. area)', fontsize = 18)
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(7.0, 5.3)
#ax.text(0.07, 0.77, leg, fontsize=13, color = 'k', transform=ax.transAxes)
adjustFigAspect(fig,aspect=1.0)
plt.gcf().subplots_adjust(bottom=0.19, left = 0.15)
plt.ylim(0.4,1.8)
plt.xlim(0, radius + 2)
plt.yticks(np.arange(0.5, 2.0, 0.5))
plt.xticks(np.arange(10, radius + 1, 10))
plt.savefig(savePath + "%s_dist_vs_aA_r_%s_z_%s_%s_fixed_range.pdf"%(label, radius, zMin, region), format='pdf', dpi=300)
plt.close()


colours = (np.array(all_dist) - np.min(all_dist))/(np.max(all_dist) - np.min(all_dist))
m, b, r, p, std_err = linregress(all_times, np.array(all_vols)/meanV)
print m, b, r, p
leg = "$N$ = " + str(len(all_vols)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
fig =plt.figure()   
ax = fig.add_subplot(1,1,1)
ax.scatter(all_times, np.array(all_vols)/meanV, c= colours, alpha = 0.3, s=50,  marker ="o", edgecolors='None')
#plt.(0.0, 4.0)
#ax.plot(all_dist, m * np.array(all_dist) + b, 'k-.')
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(5, 3.5)
ax.text(0.07, 0.8, leg, fontsize=6, color = 'k', transform=ax.transAxes)
adjustFigAspect(fig,aspect=1.0)
plt.gcf().subplots_adjust(bottom=0.19, left = 0.19)
plt.xlabel('$t$  (h)')
plt.ylabel('$V/\mu_V$')
plt.ylim(-0.1,2.5)
plt.xlim(t_begin, times[-1] + 2)
plt.yticks(np.arange(-.0, 3, 0.5))
plt.xticks(np.arange(t_begin + 10, times[-1] + 1, 10))
plt.savefig(savePath + "%s_t_vs_V_r_%s_z_%s_tb_%s_%s_fixed_range.pdf"%(label, radius, zMin, t_begin, region), format='pdf', dpi=300)
plt.close()


sk = skew(np.array(all_vols)/meanV); skT = skewtest(np.array(all_vols)/meanV)
leg = "$N$ = " + str(len(all_vols)) + "\n" + "skewness, $p$ = " + str(round(sk,2)) + ', ' + str(round_sig(skT[1],2)) +  "\n"
fig =plt.figure() 
ax = fig.add_subplot(1,1,1)  
bins = [0,0.2,0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2,2.2]
plt.hist(np.array(all_vols)/meanV, bins, alpha=0.5, color = 'c')
plt.axvline(np.median(np.array(all_vols)/meanV), color='r', linestyle='dashed', linewidth=2)
plt.xlim(0, 2.2)
plt.xlabel('$V/\mu_V$')
#plt.ylabel('frequency')
#adjustFigAspect(fig,aspect=0.8)
ax.text(0.07, 0.85, leg, fontsize=6, color = 'k', transform=ax.transAxes)
fig.set_size_inches(5.5, 4.0)
plt.gcf().subplots_adjust(bottom=0.12, left = 0.12)
plt.yticks(np.arange(0.0, 401, 100))
plt.savefig(savePath + "%s_hist_V_r_%s_z_%s_tb_%s_%s_fixed_range.pdf"%(label, radius, zMin, t_begin, region), format='pdf', dpi=300)
plt.close()


sk = skew(np.array(all_opA)/meanopA); skT = skewtest(np.array(all_opA)/meanopA)
leg = "$N$ = " + str(len(all_opA)) + "\n" + "skew., $p$ = " + str(round(sk,2)) + ', ' + str(round_sig(skT[1],2)) +  "\n"
fig =plt.figure() 
ax = fig.add_subplot(1,1,1)  
bins = [0,0.2,0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2,2.2]
plt.hist(np.array(all_opA)/meanopA, bins, alpha=0.5, color = 'c')
plt.axvline(np.median(np.array(all_opA)/meanopA), color='r', linestyle='dashed', linewidth=2)
plt.xlim(0, 2.2)
plt.xlabel('outer peri. area/mean(outer peri. area)', fontsize =20)
#plt.ylabel('frequency')
#adjustFigAspect(fig,aspect=0.8)
ax.text(0.07, 0.85, leg, fontsize=13, color = 'k', transform=ax.transAxes)
fig.set_size_inches(6.5, 5.2)
plt.gcf().subplots_adjust(bottom=0.15, left = 0.15)
plt.yticks(np.arange(100, 601, 100))
plt.savefig(savePath + "%s_hist_opA_r_%s_z_%s_tb_%s_%s_fixed_range.pdf"%(label, radius, zMin, t_begin, region), format='pdf', dpi=300)
plt.close()


sk = skew(np.array(all_num_L1neighs)); skT = skewtest(np.array(all_num_L1neighs))
leg = "$N$ = " + str(len(all_num_L1neighs)) + "\n" + "skew., $p$ = " + str(round(sk,2)) + ', ' + str(round_sig(skT[1],2)) +  "\n"
#leg = "n = " + str(len(all_num_L1neighs))
fig =plt.figure() 
ax = fig.add_subplot(1,1,1)  
bins = [2.5,3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5]
plt.hist(np.array(all_num_L1neighs), bins, alpha=0.5, color = 'c')
plt.axvline(np.median(np.array(all_num_L1neighs)), color='r', linestyle='dashed', linewidth=2)
plt.xlim(3.5, 10)
plt.xlabel('$N_{neigh}$', fontsize = 20)
#plt.ylabel('frequency')
#adjustFigAspect(fig,aspect=0.8)
ax.text(0.52, 0.8, leg, fontsize=13, color = 'k', transform=ax.transAxes)
fig.set_size_inches(6.5, 5.2)
plt.yticks(np.arange(0.0, 1001, 200))
plt.gcf().subplots_adjust(bottom=0.15, left = 0.15)
plt.savefig(savePath + "%s_hist_Nneighs_r_%s_z_%s_tb_%s_%s_fixed_range.pdf"%(label, radius, zMin, t_begin, region), format='pdf', dpi=300)
plt.close()


Asplit = [ np.amin(all_num_L1neighs), 5,6, 7, 8, np.amax(all_num_L1neighs)]
Abins = [[] for k in range(len(Asplit)-1)]; dABins = [[] for k in range(len(Asplit)-1)]
for k in range(len(all_num_L1neighs)):
    for l in range(len(Asplit)-1):
        if (all_num_L1neighs[k] >= Asplit[l] and all_num_L1neighs[k] < Asplit[l + 1]): 
            Abins[l].append(all_num_L1neighs[k]); dABins[l].append(all_vols[k]/meanV)   
x = [np.median(Abins[k]) for k in range(len(Abins))]; y = [np.median(dABins[k]) for k in range(len(dABins))]
yerr = [[y[k] - np.percentile(dABins[k], 25) for k in range(len(dABins))], [np.percentile(dABins[k], 75) - y[k] for k in range(len(dABins))]]           

m, b, r, p, std_err = linregress(x, y)
print m, b, r, p
leg = "$N$ = " + str(len(all_vols)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
fig =plt.figure()   
ax = fig.add_subplot(1,1,1)
ax.scatter(all_num_L1neighs, np.array(all_vols)/meanV, color ="0.55", alpha = 0.3,  marker ="o", edgecolors='None')
#plt.(0.0, 4.0)
ax.plot(all_num_L1neighs, m * np.array(all_num_L1neighs) + b, 'b-.')
ax.errorbar(x, y, yerr =  yerr, linestyle='-', color = 'k', marker='o')
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(5, 3.5)
ax.text(0.07, 0.8, leg, fontsize=6, color = 'k', transform=ax.transAxes)
adjustFigAspect(fig,aspect=1.0)
plt.gcf().subplots_adjust(bottom=0.15, left = 0.15)
plt.xlabel('$N_{neigh}$')
plt.ylabel('$V/\mu_V$')
plt.ylim(-0.1,2.5)
plt.xlim(2, 11)
#plt.yticks(np.arange(-.0, 3, 0.5))
#plt.xticks(np.arange(5, radius + 1, 5))
plt.savefig(savePath + "Nneigh_vs_V_r_%s_z_%s_%s_fixed_range.pdf"%(radius, zMin, region), format='pdf', dpi=300)
plt.close()


dABins = [[] for k in range(len(Asplit)-1)]
for k in range(len(all_num_L1neighs)):
    for l in range(len(Asplit)-1):
        if (all_num_L1neighs[k] >= Asplit[l] and all_num_L1neighs[k] < Asplit[l + 1]): 
             dABins[l].append(all_opA[k]/meanopA)   
x = [np.median(Abins[k]) for k in range(len(Abins))]; y = [np.median(dABins[k]) for k in range(len(dABins))]
yerr = [[y[k] - np.percentile(dABins[k], 25) for k in range(len(dABins))], [np.percentile(dABins[k], 75) - y[k] for k in range(len(dABins))]]           

m, b, r, p, std_err = linregress(all_num_L1neighs, np.array(all_opA)/meanopA)
print m, b, r, p
leg = "$N$ = " + str(len(all_opA)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
fig =plt.figure()   
ax = fig.add_subplot(1,1,1)
#ax.scatter(all_num_L1neighs, np.array(all_opA)/meanopA, color ="0.55", alpha = 0.1,  marker ="o", edgecolors='None')
ax.scatter(all_num_L1neighs, np.array(all_opA)/meanopA, color = 'c', s = 120, alpha = 0.02,  marker ="o", edgecolors = 'None')
#plt.(0.0, 4.0)
ax.plot(np.arange(4, 9, 1), m * np.arange(4, 9, 1) +b, 'r-')
ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(4.2, 4.2)
#ax.text(0.07, 0.8, leg, fontsize=6, color = 'k', transform=ax.transAxes)
adjustFigAspect(fig,aspect=1.0)
plt.xlabel('$N_{neigh}$', fontsize = 18)
plt.ylabel('$A_{op}/\mu_{A_{op}}$', fontsize = 18)
plt.ylim(0.3,1.7)
plt.xlim(3.0, 9.0)
plt.yticks(np.arange(0.5, 2.0, 0.5))
#plt.xticks(np.arange(5, radius + 1, 5))
plt.gcf().subplots_adjust(bottom=0.18, left = 0.22)
plt.savefig(savePath + "Nneigh_vs_Aop_r_%s_z_%s_%s_fixed_range.pdf"%(radius, zMin, region), format='pdf', dpi=300)
plt.close()


fig = plt.figure()
ax = fig.add_subplot(111)
sm.graphics.violinplot(vols_t, positions = timesLabel, labels = timesLabel2, ax=ax, show_boxplot = True,
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs', 'violin_lw':0,
                                  'label_fontsize':'large',
                                  'label_rotation':30, 'violin_width':3.0})
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = ( 0 - times[0], 16 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = ( 24 - times[0], 24 + 16 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = ( 2*24 - times[0], 2*24 + 16 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = ( 3*24 - times[0], 3*24 + 16 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest', alpha =0.1, extent = (16 - times[0], 24 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest', alpha =0.1, extent = (16 + 24 - times[0], 24 + 24 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest',alpha =0.1, extent = (16 + 2*24 - times[0], 24 + 2 *24 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest',alpha =0.1, extent = (16 + 3*24 - times[0], 24 +  3*24 - times[0], 0, 2*350))
ax.plot(np.arange(0, times[-1] -times[0] + 5, 1), meanV + np.zeros(len(np.arange(0, times[-1] -times[0] + 5, 1))),'b--')
plt.ylabel('vol. ($\mu$m$^3$)', fontsize = 20)
plt.xlim(-5, times[-1] - times[0] + 5)
plt.ylim(0,350)
plt.xlabel('time (h)', fontsize = 20)
#adjustFigAspect(fig,aspect=0.8)
ax.set_aspect('auto')
plt.gcf().subplots_adjust(bottom=0.22, left = 0.15)
fig.set_size_inches(9.0, 4.3)
plt.savefig(savePath + "%s_t_vs_v_violin_box.pdf"%(label), format='pdf', dpi=300)
plt.close()


fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on = False)
sm.graphics.violinplot(vols_t, positions = np.array(times) - times[0], labels = timesLabel2, ax=ax, show_boxplot = False,
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs', 'violin_lw':0,
                                  'label_fontsize':'large',
                                  'label_rotation':0, 'violin_width':2.6, 'violin_alpha':0.7})
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = ( 0 - times[0], 16 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = ( 24 - times[0], 24 + 16 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = ( 2*24 - times[0], 2*24 + 16 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = ( 3*24 - times[0], 3*24 + 16 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest', alpha =0.1, extent = (16 - times[0], 24 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest', alpha =0.1, extent = (16 + 24 - times[0], 24 + 24 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest',alpha =0.1, extent = (16 + 2*24 - times[0], 24 + 2 *24 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest',alpha =0.1, extent = (16 + 3*24 - times[0], 24 +  3*24 - times[0], 0, 2*350))
ax.plot(np.arange(0, times[-1] -times[0] + 5, 1), meanV + np.zeros(len(np.arange(0, times[-1] -times[0] + 5, 1))),'b--')
plt.ylabel('$V$ ($\mu$m$^3$)', fontsize = 25)
plt.xlim(-5, times[-1] -times[0] + 5)
plt.ylim(25,300)
plt.xlabel('$t$ (h)', fontsize = 25)
#adjustFigAspect(fig,aspect=0.8)
ax.set_aspect('auto')
plt.gcf().subplots_adjust(bottom=0.17, left = 0.1)
fig.set_size_inches(11.0, 4.2)
plt.savefig(savePath + "%s_t_vs_v_violin.pdf"%label, format='pdf', dpi=300)
plt.close()


fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = ( 0 - times[0], 16 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = ( 24 - times[0], 24 + 16 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = ( 2*24 - times[0], 2*24 + 16 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = ( 3*24 - times[0], 3*24 + 16 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest', alpha =0.1, extent = (16 - times[0], 24 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest', alpha =0.1, extent = (16 + 24 - times[0], 24 + 24 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest',alpha =0.1, extent = (16 + 2*24 - times[0], 24 + 2 *24 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest',alpha =0.1, extent = (16 + 3*24 - times[0], 24 +  3*24 - times[0], 0, 2*350))
ax.plot(np.array(times[1:]) - times[0], np.array(numberDivisions)/(np.array(times[1:]) - np.array(times[:-1])),'b--')
plt.ylabel('division rate in central zone (h$^{-1}$)')
plt.xlim(times[0]-5, times[-1]-times[0] + 5)
plt.ylim(0, 19)
plt.xlabel('$t$ (h)')
#adjustFigAspect(fig,aspect=0.8)
ax.set_aspect('auto')
plt.gcf().subplots_adjust(bottom=0.14, left = 0.1)
fig.set_size_inches(6.5, 3.5)
plt.savefig(savePath + "%s_t_vs_num_divisions_%s.pdf"%(label,radius), format='pdf', dpi=300)
plt.close()

print np.array(numberDivisionsPerCell)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = ( 0 - times[0], 16 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = ( 24 - times[0], 24 + 16 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = ( 2*24 - times[0], 2*24 + 16 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.spring, interpolation = 'nearest', alpha =0.1, extent = ( 3*24 - times[0], 3*24 + 16 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest', alpha =0.1, extent = (16 - times[0], 24 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest', alpha =0.1, extent = (16 + 24 - times[0], 24 + 24 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest',alpha =0.1, extent = (16 + 2*24 - times[0], 24 + 2 *24 - times[0], 0, 2*350))
ax.imshow([[50,50], [100,100]], cmap=plt.cm.Blues, interpolation = 'nearest',alpha =0.1, extent = (16 + 3*24 - times[0], 24 +  3*24 - times[0], 0, 2*350))
ax.plot(np.array(times[1:]) - times[0], np.array(numberDivisionsPerCell)/(np.array(times[1:]) - np.array(times[:-1])),'b--')
plt.ylabel('division rate per cell (h$^{-1}$)', fontsize = 20)
plt.xlim(times[0]-5, times[-1]-times[0] + 5)
plt.ylim(0, 0.14)
plt.xlabel( 'time (h)', fontsize = 20)
#adjustFigAspect(fig,aspect=0.8)
ax.set_aspect('auto')
plt.gcf().subplots_adjust(bottom=0.16, left = 0.16)
fig.set_size_inches(9, 4.5)
plt.savefig(savePath + "%s_t_vs_num_divisions_per_cell_%s.pdf"%(label,radius), format='pdf', dpi=300)
plt.close()
