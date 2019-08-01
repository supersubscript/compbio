import cPickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.stats import pearsonr, skew, skewtest
import statsmodels.api as sm

radius = '10'
zThres = '8'
region = 'L1'

plt.rcParams['xtick.labelsize'] = 13; plt.rcParams['ytick.labelsize'] = 13

#volumes over cell cycle with last volume being the total volume of two daughters upon division
filePathsVols = [
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant2/statistics/birth_div_vols_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant4/statistics/birth_div_vols_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant13/statistics/birth_div_vols_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant15/statistics/birth_div_vols_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant18/statistics/birth_div_vols_r_%s_z_%s_%s.out"%(radius, zThres, region)
]

filePathsNormVols = [
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant2/statistics/birth_div_norm_vols_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant4/statistics/birth_div_norm_vols_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant13/statistics/birth_div_norm_vols_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant15/statistics/birth_div_norm_vols_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant18/statistics/birth_div_norm_vols_r_%s_z_%s_%s.out"%(radius, zThres, region)
]

filePathsNucVols = [
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant2/statistics/birth_div_nuc_vols_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant4/statistics/birth_div_nuc_vols_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant13/statistics/birth_div_nuc_vols_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant15/statistics/birth_div_nuc_vols_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant18/statistics/birth_div_nuc_vols_r_%s_z_%s_%s.out"%(radius, zThres, region)
]

filePathsNormNucVols = [
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant2/statistics/birth_div_norm_nuc_vols_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant4/statistics/birth_div_norm_nuc_vols_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant13/statistics/birth_div_norm_nuc_vols_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant15/statistics/birth_div_norm_nuc_vols_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant18/statistics/birth_div_norm_nuc_vols_r_%s_z_%s_%s.out"%(radius, zThres, region)
]

filePathsTimes = [
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant2/statistics/birth_div_times_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant4/statistics/birth_div_times_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant13/statistics/birth_div_times_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant15/statistics/birth_div_times_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant18/statistics/birth_div_times_r_%s_z_%s_%s.out"%(radius, zThres, region)
]

filePathsDists = [
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant2/statistics/birth_div_dist_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant4/statistics/birth_div_dist_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant13/statistics/birth_div_dist_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant15/statistics/birth_div_dist_r_%s_z_%s_%s.out"%(radius, zThres, region),
"/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant18/statistics/birth_div_dist_r_%s_z_%s_%s.out"%(radius, zThres, region)
]

def round_sig(x, sig=2):
    if x ==0 : return 0.0
    else: return np.round(x, sig-int(np.floor(np.log10(x)))-1)

def create_bins(n, arr_x, arr_y):
    Asplit = [np.percentile(arr_x, float(i)/n*100) for i in range(n + 1)]
    Abins = [[] for k in range(len(Asplit)-1)]; dABins = [[] for k in range(len(Asplit)-1)]
    for k in range(len(arr_x)):
        for l in range(len(Asplit)-1):
            if (arr_x[k] >= Asplit[l] and arr_x[k] < Asplit[l + 1]): 
                Abins[l].append(arr_x[k]); dABins[l].append(arr_y[k]) 
      
    x = [np.median(Abins[k]) for k in range(len(Abins))]; y = [np.median(dABins[k]) for k in range(len(dABins))]
    yerr = [[y[k] - np.percentile(dABins[k], 25) for k in range(len(dABins))], [np.percentile(dABins[k], 75) - y[k] for k in range(len(dABins))]]
    return x, y, yerr     

all_vols = []
all_nuc_vols = []
all_norm_vols = []
all_norm_nuc_vols = []
all_division_flags = []
all_dist = []
all_times = []
all_opAs = []
all_plants = []

for fl in filePathsNucVols:
    with open(fl, "r") as infile:
        nuc_volumes_in = [np.array(map(str, line.split(','))) for line in infile]
        for vol in nuc_volumes_in:
            for i in range(len(vol) -1):
                vol[i] = vol[i][1:]
            vol[-1] = vol[-1][1:-2]
    nuc_volumes = [map(float, arr) for arr in nuc_volumes_in]
    all_nuc_vols.append(nuc_volumes)

for fl in filePathsNormNucVols:
    with open(fl, "r") as infile:
        nuc_volumes_in = [np.array(map(str, line.split(','))) for line in infile]
        for vol in nuc_volumes_in:
            for i in range(len(vol) -1):
                vol[i] = vol[i][1:]
            vol[-1] = vol[-1][1:-2]
    nuc_volumes = [map(float, arr) for arr in nuc_volumes_in]
    all_norm_nuc_vols.append(nuc_volumes)

for fl in filePathsVols:
    with open(fl, "r") as infile:
        volumes_in = [np.array(map(str, line.split(','))) for line in infile]
        for vol in volumes_in:
            for i in range(len(vol) -1):
                vol[i] = vol[i][1:]
            vol[-1] = vol[-1][1:-2]
    volumes = [map(float, arr) for arr in volumes_in]
    all_vols.append(volumes)

for fl in filePathsNormVols:
    with open(fl, "r") as infile:
        volumes_in = [np.array(map(str, line.split(','))) for line in infile]
        for vol in volumes_in:
            for i in range(len(vol) -1):
                vol[i] = vol[i][1:]
            vol[-1] = vol[-1][1:-2]
    volumes = [map(float, arr) for arr in volumes_in]
    all_norm_vols.append(volumes)

for fl in filePathsTimes:
    with open(fl, "r") as infile:
        times_in = [np.array(map(str, line.split(','))) for line in infile]
        for tm in times_in:
            for i in range(len(tm) -1):
                tm[i] = tm[i][1:]
            tm[-1] = tm[-1][1:-2]
    times = [map(float, arr) for arr in times_in]
    all_times.append(times)

num_plants = 0
for fl in filePathsDists:
    with open(fl, "r") as infile:
        dist_in = [np.array(map(str, line.split(','))) for line in infile]
        for ds in dist_in:
            for i in range(len(ds) -1):
                ds[i] = ds[i][1:]
            ds[-1] = ds[-1][1:-2]
    dist = [map(float, arr) for arr in dist_in]
    all_plants.append([])
    for entry in dist:
        all_plants[-1].append(num_plants)
    all_dist.append(dist)
    num_plants = num_plants + 1 

print all_dist
print all_plants
    
#normalize by means
birth_nuc_vols = [cc_vol[0] for vol_plant in all_nuc_vols for cc_vol in vol_plant]
div_nuc_vols =  [cc_vol[-1] for vol_plant in all_nuc_vols for cc_vol in vol_plant]
av_dist =  [np.mean(all_dist[j][i]) for j in range(len(all_dist)) for i in range(len(all_dist[j]))]

birth_norm_nuc_vols = [cc_vol[0] for vol_plant in all_norm_nuc_vols for cc_vol in vol_plant if cc_vol[0]>0.05 and cc_vol[0] <0.45 and len(cc_vol) > 1]
div_norm_nuc_vols =  [cc_vol[-1] for vol_plant in all_norm_nuc_vols for cc_vol in vol_plant if cc_vol[0]>0.05 and cc_vol[0] <0.45 and len(cc_vol) > 1]
b4_div_norm_nuc_vols =  [cc_vol[-2] for vol_plant in all_norm_nuc_vols for cc_vol in vol_plant if cc_vol[0]>0.05 and cc_vol[0] <0.45 and len(cc_vol) > 1]
av_norm_dist =  [np.mean(all_dist[j][i]) for j in range(len(all_dist)) for i in range(len(all_dist[j])) if all_norm_nuc_vols[j][i][0]>0.05 and all_norm_nuc_vols[j][i][0] <0.45 and len(all_norm_nuc_vols[j][i]) > 1]

birth_vols = [cc_vol[0] for vol_plant in all_vols for cc_vol in vol_plant]
div_vols =  [cc_vol[-1] for vol_plant in all_vols for cc_vol in vol_plant]
av_dist_vol =  [np.mean(all_dist[j][i]) for j in range(len(all_dist)) for i in range(len(all_dist[j]))]
int_div_times_vol = [all_times[j][i][-1] - all_times[j][i][0] for j in range(len(all_times)) for i in range(len(all_times[j]))]


cc_phase = [(all_times[j][i][k] - all_times[j][i][0])/(all_times[j][i][-1] - all_times[j][i][0]) for j in range(len(all_times)) for i in range(len(all_times[j])) for k in range(len(all_times[j][i]))]
print "len: cc_phase", len(cc_phase)

nuc_vol_cc_phase = [cc_vol[k] for vol_plant in all_nuc_vols for cc_vol in vol_plant for k in range(len(cc_vol))]
print "len: nuc_vol_cc_phase", len(nuc_vol_cc_phase)

vol_cc_phase = [cc_vol[k] for vol_plant in all_vols for cc_vol in vol_plant for k in range(len(cc_vol))]
print "len: nuc_vol_cc_phase", len(nuc_vol_cc_phase)



birth_nuc_vols_no_out = [all_nuc_vols[j][i][0] for j in range(len(all_nuc_vols)) for i in range(len(all_nuc_vols[j])) if all_nuc_vols[j][i][0] > 10 and all_nuc_vols[j][i][0] < 60 and all_nuc_vols[j][i][-2] >  10 and len(all_nuc_vols[j][i]) > 1] 
div_nuc_vols_no_out =  [all_nuc_vols[j][i][-2] for j in range(len(all_nuc_vols)) for i in range(len(all_nuc_vols[j])) if all_nuc_vols[j][i][0] > 10 and all_nuc_vols[j][i][0] < 60 and all_nuc_vols[j][i][-2] >  10 and len(all_nuc_vols[j][i]) > 1] 
av_dist_no_out =  [np.mean(all_dist[j][i]) for j in range(len(all_dist)) for i in range(len(all_dist[j])) if all_nuc_vols[j][i][0] > 10 and all_nuc_vols[j][i][0] < 60 and all_nuc_vols[j][i][-2] >  10 and len(all_nuc_vols[j][i]) > 1]
int_div_times_no_out = [all_times[j][i][-1] - all_times[j][i][0] for j in range(len(all_times)) for i in range(len(all_times[j])) if all_nuc_vols[j][i][0] > 10 and all_nuc_vols[j][i][0] < 60 and all_nuc_vols[j][i][-2] >  10 and len(all_nuc_vols[j][i]) > 1]
all_plants_no_out = [all_plants[j][i] for j in range(len(all_plants)) for i in range(len(all_plants[j])) if all_nuc_vols[j][i][0] > 10 and all_nuc_vols[j][i][0] < 60 and all_nuc_vols[j][i][-2] >  10 and len(all_nuc_vols[j][i]) > 1] 

birth_norm_cyto =  [all_norm_vols[j][i][0] - all_norm_nuc_vols[j][i][0] for j in range(len(all_norm_vols)) for i in range(len(all_norm_vols[j])) if all_norm_nuc_vols[j][i][0]>0.05 and all_norm_nuc_vols[j][i][0] <0.45 and len(all_norm_nuc_vols[j][i]) > 1]
div_norm_cyto =  [all_norm_vols[j][i][-1] - all_norm_nuc_vols[j][i][-1] for j in range(len(all_norm_vols)) for i in range(len(all_norm_vols[j])) if all_norm_nuc_vols[j][i][0]>0.05 and all_norm_nuc_vols[j][i][0] <0.45 and len(all_norm_nuc_vols[j][i]) > 1]

int_div_times = [all_times[j][i][-1] - all_times[j][i][0] for j in range(len(all_times)) for i in range(len(all_times[j])) if all_norm_nuc_vols[j][i][0]>0.05 and all_norm_nuc_vols[j][i][0] <0.45 and len(all_norm_nuc_vols[j][i]) > 1]


#N######UCLEAR VOLUMES
x, y, yerr = create_bins(2, birth_norm_nuc_vols, div_norm_nuc_vols);
fig= plt.figure()  
ax = fig.add_subplot(1, 1, 1)
s = ax.scatter(birth_norm_nuc_vols, div_norm_nuc_vols, s = 100, c=av_norm_dist, marker ="o", edgecolors='None', alpha = 0.3) 
s.set_clim([np.amin(av_norm_dist), np.amax(av_norm_dist)])
cbar = fig.colorbar(s, ticks = [2,4,6,8], shrink = 0.6, pad = 0.03)
cbar.ax.set_ylabel('distance from $O$ ($\mu$m)', fontsize = 10)
cbar.ax.set_aspect(10)
#ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
m, b = np.polyfit(birth_norm_nuc_vols, div_norm_nuc_vols, 1); r, p = pearsonr(birth_norm_nuc_vols, div_norm_nuc_vols)
ax.plot(birth_norm_nuc_vols, 1.6 * np.array(birth_norm_nuc_vols) + 0, linestyle = '-', color = 'r', linewidth = 1.3 )
ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
leg = "$N$ = " + str(len(birth_norm_nuc_vols)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
ax.text(0.54, 0.79, leg, fontsize=11, color = 'k', transform=ax.transAxes)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)
fig.set_size_inches(6.5, 6.0)
plt.gcf().subplots_adjust(bottom=0.15, left = 0.13)
plt.savefig("/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_nuclear_data/BirthNormNucVolVsDivNormNucVol_r_%s_z_%s_%s.pdf"%(radius, zThres,region), format='pdf', dpi=200)
plt.close(fig)


x, y, yerr = create_bins(2,birth_nuc_vols, div_nuc_vols);
fig= plt.figure()  
ax = fig.add_subplot(1, 1, 1)
s = ax.scatter(birth_nuc_vols, div_nuc_vols, s = 100, c=av_dist, marker ="o", edgecolors='None', alpha = 0.3) 
s.set_clim([np.amin(av_dist), np.amax(av_dist)])
cbar = fig.colorbar(s, ticks = [2,4,6,8], shrink = 0.6, pad = 0.03)
cbar.ax.set_ylabel('distance from $O$ ($\mu$m)', fontsize = 10)
cbar.ax.set_aspect(10)
#ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
m, b = np.polyfit(birth_nuc_vols, div_nuc_vols, 1); r, p = pearsonr(birth_nuc_vols, div_nuc_vols)
ax.plot(birth_nuc_vols, 1.8 * np.array(birth_nuc_vols) + 0, linestyle = '-', color = 'r', linewidth = 1.3 )
ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
leg = "$N$ = " + str(len(birth_nuc_vols)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
ax.text(0.54, 0.79, leg, fontsize=11, color = 'k', transform=ax.transAxes)
ax.set_xlabel(' birth nuclear vol. ($\mu$m$^3$)', fontsize =20)
ax.set_ylabel(' division nuclear vol. ($\mu$m$^3$)', fontsize = 20)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)
fig.set_size_inches(6.5, 6.0)
plt.gcf().subplots_adjust(bottom=0.15, left = 0.15)
plt.savefig("/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_nuclear_data/BirthNucVolVsDivNucVol_r_%s_z_%s_%s.pdf"%(radius, zThres,region), format='pdf', dpi=200)
plt.close(fig)


x, y, yerr = create_bins(2, birth_nuc_vols_no_out, div_nuc_vols_no_out);
fig= plt.figure()  
ax = fig.add_subplot(1, 1, 1)
s = ax.scatter(birth_nuc_vols_no_out, div_nuc_vols_no_out, s = 100, c=av_dist_no_out, marker ="o", edgecolors='None', alpha = 0.3) 
s.set_clim([np.amin(av_dist_no_out), np.amax(av_dist_no_out)])
cbar = fig.colorbar(s, ticks = [2,4,6,8], shrink = 0.6, pad = 0.03)
cbar.ax.set_ylabel('distance from $O$ ($\mu$m)', fontsize = 10)
cbar.ax.set_aspect(10)
#ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
m, b = np.polyfit(birth_nuc_vols_no_out, div_nuc_vols_no_out, 1); r, p = pearsonr(birth_nuc_vols_no_out, div_nuc_vols_no_out)
ax.plot(birth_nuc_vols_no_out, 1.8 * np.array(birth_nuc_vols_no_out) + 0, linestyle = '-', color = 'r', linewidth = 1.3 )
ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
leg = "$N$ = " + str(len(birth_nuc_vols_no_out)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
ax.text(0.54, 0.79, leg, fontsize=11, color = 'k', transform=ax.transAxes)
ax.set_xlabel(' birth nuclear vol. ($\mu$m$^3$)', fontsize =20)
ax.set_ylabel(' division nuclear vol. ($\mu$m$^3$)', fontsize = 20)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)
fig.set_size_inches(6.5, 6.0)
plt.gcf().subplots_adjust(bottom=0.15, left = 0.15)
plt.savefig("/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_nuclear_data/BirthNucVolVsDivNucVol_noOutliers_r_%s_z_%s_%s_colDist.pdf"%(radius, zThres,region), format='pdf', dpi=200)
plt.close(fig)

fig= plt.figure()  
ax = fig.add_subplot(1, 1, 1)
ax.scatter(birth_nuc_vols_no_out, div_nuc_vols_no_out, s = 100, c=all_plants_no_out, marker ="o", edgecolors='None', alpha = 0.3) 
#ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
m, b = np.polyfit(birth_nuc_vols_no_out, div_nuc_vols_no_out, 1); r, p = pearsonr(birth_nuc_vols_no_out, div_nuc_vols_no_out)
ax.plot(birth_nuc_vols_no_out, 1.8 * np.array(birth_nuc_vols_no_out) + 0, linestyle = '-', color = 'r', linewidth = 1.3 )
ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
leg = "$N$ = " + str(len(birth_nuc_vols_no_out)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
ax.text(0.54, 0.79, leg, fontsize=11, color = 'k', transform=ax.transAxes)
ax.set_xlabel(' birth nuclear vol. ($\mu$m$^3$)', fontsize =20)
ax.set_ylabel(' division nuclear vol. ($\mu$m$^3$)', fontsize = 20)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)
fig.set_size_inches(6.5, 6.0)
plt.gcf().subplots_adjust(bottom=0.15, left = 0.15)
plt.savefig("/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_nuclear_data/BirthNucVolVsDivNucVol_noOutliers_r_%s_z_%s_%s_colPlant.pdf"%(radius, zThres,region), format='pdf', dpi=200)
plt.close(fig)

x, y, yerr = create_bins(2,birth_nuc_vols_no_out, int_div_times_no_out);
fig= plt.figure()  
ax = fig.add_subplot(1, 1, 1)
s = ax.scatter(birth_nuc_vols_no_out, int_div_times_no_out, s = 100, c=av_dist_no_out, marker ="o", edgecolors='None', alpha = 0.3) 
s.set_clim([np.amin(av_dist_no_out), np.amax(av_dist_no_out)])
cbar = fig.colorbar(s, ticks = [2,4,6,8], shrink = 0.6, pad = 0.03)
cbar.ax.set_ylabel('distance from $O$ ($\mu$m)', fontsize = 10)
cbar.ax.set_aspect(10)
#ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
m, b = np.polyfit(birth_nuc_vols_no_out, int_div_times_no_out, 1); r, p = pearsonr(birth_nuc_vols_no_out, int_div_times_no_out)
#ax.plot(birth_nuc_vols_no_out, 1.8 * np.array(birth_nuc_vols_no_out) + 0, linestyle = '-', color = 'r', linewidth = 1.3 )
ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
leg = "$N$ = " + str(len(birth_nuc_vols_no_out)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
ax.text(0.54, 0.79, leg, fontsize=11, color = 'k', transform=ax.transAxes)
ax.set_xlabel(' birth nuclear vol. ($\mu$m$^3$)', fontsize =20)
ax.set_ylabel(' inter-division time', fontsize = 20)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)
fig.set_size_inches(6.5, 6.0)
plt.gcf().subplots_adjust(bottom=0.15, left = 0.15)
plt.savefig("/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_nuclear_data/BirthNucVolVsIntDivTimes_noOutliers_r_%s_z_%s_%s.pdf"%(radius, zThres,region), format='pdf', dpi=200)
plt.close(fig)

x, y, yerr = create_bins(2,birth_norm_nuc_vols, b4_div_norm_nuc_vols);
fig= plt.figure()  
ax = fig.add_subplot(1, 1, 1)
ax.scatter(birth_norm_nuc_vols, b4_div_norm_nuc_vols, s = 100, color='c', marker ="o", edgecolors='None', alpha = 0.3) 
#ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
m, b = np.polyfit(birth_norm_nuc_vols, b4_div_norm_nuc_vols, 1); r, p = pearsonr(birth_norm_nuc_vols, b4_div_norm_nuc_vols)
ax.plot(birth_norm_nuc_vols, 1.5 * np.array(birth_norm_nuc_vols) + 0, linestyle = '-', color = 'r', linewidth = 1.3 )
ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
leg = "$N$ = " + str(len(birth_norm_nuc_vols)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
ax.text(0.54, 0.79, leg, fontsize=11, color = 'k', transform=ax.transAxes)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)
fig.set_size_inches(6.5, 6.0)
plt.gcf().subplots_adjust(bottom=0.15, left = 0.13)
plt.savefig("/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_nuclear_data/BirthNormNucVolVsB4DivNormNucVol_r_%s_z_%s_%s.pdf"%(radius, zThres,region), format='pdf', dpi=200)
plt.close(fig)

x, y, yerr = create_bins(2, birth_norm_nuc_vols, np.array(div_norm_nuc_vols) - np.array(birth_norm_nuc_vols));
fig= plt.figure()  
ax = fig.add_subplot(1, 1, 1)
ax.scatter(birth_norm_nuc_vols, np.array(div_norm_nuc_vols) - np.array(birth_norm_nuc_vols), s = 100, color='c', marker ="o", edgecolors='None', alpha = 0.3) 
#ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
m, b = np.polyfit(birth_norm_nuc_vols, np.array(div_norm_nuc_vols) - np.array(birth_norm_nuc_vols), 1); r, p = pearsonr(birth_norm_nuc_vols, np.array(div_norm_nuc_vols) - np.array(birth_norm_nuc_vols))
ax.plot(birth_norm_nuc_vols, m * np.array(birth_norm_nuc_vols) + b, linestyle = '-', color = 'r', linewidth = 1.3 )
ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
leg = "$N$ = " + str(len(birth_norm_nuc_vols)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
ax.text(0.54, 0.79, leg, fontsize=11, color = 'k', transform=ax.transAxes)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)
fig.set_size_inches(6.5, 6.0)
plt.gcf().subplots_adjust(bottom=0.15, left = 0.13)
plt.savefig("/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_nuclear_data/BirthNormNucVolVsIncNormNucVol_r_%s_z_%s_%s.pdf"%(radius, zThres,region), format='pdf', dpi=200)
plt.close(fig)

x, y, yerr = create_bins(2, birth_norm_nuc_vols, int_div_times)
fig= plt.figure()  
ax = fig.add_subplot(1, 1, 1)
ax.scatter(birth_norm_nuc_vols, int_div_times, s = 100, color='c', marker ="o", edgecolors='None', alpha = 0.3) 
#ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
m, b = np.polyfit(birth_norm_nuc_vols, int_div_times, 1); r, p = pearsonr(birth_norm_nuc_vols, int_div_times)
ax.plot(birth_norm_nuc_vols, m * np.array(birth_norm_nuc_vols) + b, linestyle = '-', color = 'r', linewidth = 1.3 )
ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
leg = "$N$ = " + str(len(birth_norm_nuc_vols)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
ax.text(0.54, 0.79, leg, fontsize=11, color = 'k', transform=ax.transAxes)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)
fig.set_size_inches(6.5, 6.0)
plt.gcf().subplots_adjust(bottom=0.15, left = 0.13)
plt.savefig("/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_nuclear_data/BirthNormNucVolVsIntDivTimes_r_%s_z_%s_%s.pdf"%(radius, zThres,region), format='pdf', dpi=200)
plt.close(fig)


######### CELL VOLUMES 
x, y, yerr = create_bins(4,birth_vols, div_vols);
###nuclear volunme plots
fig= plt.figure()  
ax = fig.add_subplot(1, 1, 1)
s = ax.scatter(birth_vols, div_vols, s = 100, c=av_dist_vol, marker ="o", edgecolors='None', alpha = 0.3) 
s.set_clim([np.amin(av_dist), np.amax(av_dist)])
cbar = fig.colorbar(s, ticks = [2,4,6,8], shrink = 0.6, pad = 0.03)
cbar.ax.set_ylabel('distance from $O$ ($\mu$m)', fontsize = 10)
cbar.ax.set_aspect(10)
#ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
m, b = np.polyfit(birth_vols, div_vols, 1); r, p = pearsonr(birth_vols, div_vols)
ax.plot(birth_norm_nuc_vols, 1.6 * np.array(birth_norm_nuc_vols) + 0, linestyle = '-', color = 'r', linewidth = 1.3 )
ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
leg = "$N$ = " + str(len(birth_vols)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
ax.text(0.54, 0.79, leg, fontsize=11, color = 'k', transform=ax.transAxes)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)
fig.set_size_inches(6.5, 6.0)
plt.gcf().subplots_adjust(bottom=0.15, left = 0.13)
plt.savefig("/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_nuclear_data/BirthVolVsDivVol_r_%s_z_%s_%s.pdf"%(radius, zThres,region), format='pdf', dpi=200)
plt.close(fig)


int_div_times_vol = np.array(int_div_times_vol)/np.mean(int_div_times_vol)
x, y, yerr = create_bins(2, birth_vols, int_div_times_vol)
fig= plt.figure()  
ax = fig.add_subplot(1, 1, 1)
ax.scatter(birth_vols, int_div_times_vol, s = 100, color='c', marker ="o", edgecolors='None', alpha = 0.3) 
#ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
m, b = np.polyfit(birth_vols, int_div_times_vol, 1); r, p = pearsonr(birth_vols, int_div_times_vol)
ax.plot(birth_vols, m * np.array(birth_vols) + b, linestyle = '-', color = 'r', linewidth = 1.3 )
ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
leg = "$N$ = " + str(len(birth_vols)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
ax.text(0.54, 0.79, leg, fontsize=11, color = 'k', transform=ax.transAxes)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)
fig.set_size_inches(6.5, 6.0)
plt.gcf().subplots_adjust(bottom=0.15, left = 0.13)
plt.savefig("/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_nuclear_data/BirthVolVsIntDivTimes_r_%s_z_%s_%s.pdf"%(radius, zThres,region), format='pdf', dpi=200)
plt.close(fig)

########CYTOPLASM
print "cyto"
x, y, yerr = create_bins(2,birth_norm_cyto, div_norm_cyto)
fig= plt.figure()  
ax = fig.add_subplot(1, 1, 1)
s = ax.scatter(birth_norm_cyto, div_norm_cyto, s = 100, c=av_norm_dist, marker ="o", edgecolors='None', alpha = 0.3) 
s.set_clim([np.amin(av_norm_dist), np.amax(av_norm_dist)])
cbar = fig.colorbar(s, ticks = [2,4,6,8], shrink = 0.6, pad = 0.03)
cbar.ax.set_ylabel('distance from $O$ ($\mu$m)', fontsize = 10)
cbar.ax.set_aspect(10)
#ax.set_xlim(xLimits[i]); ax.set_ylim([0.4, 2.0])
m, b = np.polyfit(birth_norm_cyto, div_norm_cyto, 1); r, p = pearsonr(birth_norm_cyto, div_norm_cyto)
ax.plot(birth_norm_cyto, m * np.array(birth_norm_cyto) + b, linestyle = '-', color = 'r', linewidth = 1.3 )
ax.errorbar(x, y, yerr =  yerr, color = 'k', fmt='o')
leg = "$N$ = " + str(len(birth_norm_cyto)) + "\n" + "sl., int. = "  + str(round(m,2)) + ", " + str(round(b,2)) + "\n" + "$R$, $p$ = "  + str(round(r,2)) + ", " + str(round_sig(p)) + "\n"
ax.text(0.54, 0.79, leg, fontsize=11, color = 'k', transform=ax.transAxes)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)
fig.set_size_inches(6.5, 6.0)
plt.gcf().subplots_adjust(bottom=0.15, left = 0.13)
plt.savefig("/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_nuclear_data/BirthNormCytoVsDivNormCyto_r_%s_z_%s_%s.pdf"%(radius, zThres,region), format='pdf', dpi=200)
plt.close(fig)

ratios = np.array([div_nuc_vols[i] for i in range(len(div_nuc_vols)) if birth_nuc_vols[i] != 0 and div_nuc_vols[i] != 0])/np.array([birth_nuc_vols[i] for i in range(len(birth_nuc_vols)) if birth_nuc_vols[i] != 0 and div_nuc_vols[i] != 0])
print "div. nuc. vol./birth nuc. vol., mean, median, cov = "
print np.mean(ratios), np.median(ratios), np.std(ratios)/np.mean(ratios)
print ratios

ratios = np.array([div_norm_nuc_vols[i] for i in range(len(div_norm_nuc_vols)) if birth_norm_nuc_vols[i] != 0 and div_norm_nuc_vols[i] != 0])/np.array([birth_norm_nuc_vols[i] for i in range(len(birth_norm_nuc_vols)) if birth_norm_nuc_vols[i] != 0 and div_norm_nuc_vols[i] != 0])
print "div. norm. nuc. vol./birth norm. nuc. vol., mean, median, cov = "
print np.mean(ratios), np.median(ratios), np.std(ratios)/np.mean(ratios)
print ratios
   
ratios = np.array([div_norm_nuc_vols[i] for i in range(len(div_norm_nuc_vols)) if b4_div_norm_nuc_vols[i] != 0 and div_norm_nuc_vols[i] != 0])/np.array([b4_div_norm_nuc_vols[i] for i in range(len(b4_div_norm_nuc_vols)) if b4_div_norm_nuc_vols[i] != 0 and div_norm_nuc_vols[i] != 0])
print "div. norm. nuc. vol./b4 div. norm. nuc. vol., mean, median, cov = "
print np.mean(ratios), np.median(ratios), np.std(ratios)/np.mean(ratios)
print ratios
 
ratios = np.array([div_norm_cyto[i] for i in range(len(div_norm_cyto)) if birth_norm_cyto[i] != 0])/np.array([birth_norm_cyto[i] for i in range(len(birth_norm_cyto)) if birth_norm_cyto[i] != 0])
print np.mean(ratios), np.median(ratios), np.std(ratios)/np.mean(ratios)
print ratios


################### GROWTH RATES

cc_frac =  [(0.5*np.array(all_times[i][j][1:-1]) + 0.5*np.array(all_times[i][j][:-2]) - all_times[i][j][0])/(all_times[i][j][-1] - all_times[i][j][0]) for i in range(len(all_times)) for j in range(len(all_times[i])) ]
cc_nuc_gr = [(np.array(all_nuc_vols[i][j][1:-1]) - np.array(all_nuc_vols[i][j][:-2]))/(0.5*np.array(all_nuc_vols[i][j][1:-1]) + 0.5*np.array(all_nuc_vols[i][j][:-2]))*1.0/(np.array(all_times[i][j][1:-1])  - np.array(all_times[i][j][:-2])) for i in range(len(all_nuc_vols)) for j in range(len(all_nuc_vols[i])) ]
cc_nuc_gr_lin = [(np.array(all_nuc_vols[i][j][1:-1]) - np.array(all_nuc_vols[i][j][:-2]))/np.mean(birth_nuc_vols)*1.0/(np.array(all_times[i][j][1:-1])  - np.array(all_times[i][j][:-2])) for i in range(len(all_nuc_vols)) for j in range(len(all_nuc_vols[i])) ]
cc_nuc_volsU = [all_nuc_vols[i][j][1:-1] for i in range(len(all_nuc_vols)) for j in range(len(all_nuc_vols[i]))]
cc_nuc_volsD = [all_nuc_vols[i][j][:-2] for i in range(len(all_nuc_vols)) for j in range(len(all_nuc_vols[i]))]
cc_vol_gr = [(np.array(all_vols[i][j][1:-1]) - np.array(all_vols[i][j][:-2]))/(0.5*np.array(all_vols[i][j][1:-1]) + 0.5*np.array(all_vols[i][j][:-2]))*1.0/(np.array(all_times[i][j][1:-1])  - np.array(all_times[i][j][:-2])) for i in range(len(all_vols)) for j in range(len(all_vols[i]))]
	
cc_frac_new = [item for sublist in cc_frac for item in sublist] 
cc_nuc_gr_new = [item for sublist in cc_nuc_gr for item in sublist] 
cc_nuc_gr_lin_new = [item for sublist in cc_nuc_gr_lin for item in sublist] 
cc_nuc_volsU = [item for sublist in cc_nuc_volsU for item in sublist] 
cc_nuc_volsD = [item for sublist in cc_nuc_volsD for item in sublist]
cc_vol_gr_new = [item for sublist in cc_vol_gr for item in sublist]

print len(cc_frac_new)
print len(cc_nuc_gr_new)
cc_frac_new_mod = [cc_frac_new[i] for i in range(len(cc_frac_new)) if cc_nuc_volsU[i] != 0 and cc_nuc_volsD[i] != 0 and not np.isinf(cc_nuc_gr_lin_new[i])] 
cc_nucgr_lin_new_mod = [cc_nuc_gr_lin_new[i] for i in range(len(cc_frac_new)) if cc_nuc_volsU[i] != 0 and cc_nuc_volsD[i] != 0 and not np.isinf(cc_nuc_gr_lin_new[i])] 

CCBins = [0, 0.25, 0.5 , 0.75, 1.0]
grBins = [[] for k in range(len(CCBins)-1)]
meanNVGR = np.mean(cc_nucgr_lin_new_mod)

for k in range(len(cc_frac_new_mod)):
    for l in range(len(CCBins)-1): 
        if (cc_frac_new_mod[k] >= CCBins[l] and cc_frac_new_mod[k] < CCBins[l + 1]): 
            grBins[l].append(cc_nucgr_lin_new_mod[k])
print grBins
labs = []
for t in range(len(CCBins)-1):
    labs.append(str(CCBins[t]) + '~' + str(CCBins[t+1]))

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
sm.graphics.violinplot(grBins, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs', 'violin_lw':0,
                                  'label_fontsize':'small',
                                  'label_rotation':30})
ax.set_xlabel(' frac. inter-division time ')
ax.set_ylabel(' nuc. vol. lin. growth rate ')
ax.set_ylim([-0.05,0.1])
fig.set_size_inches(5.0, 4.2)
plt.gcf().subplots_adjust(bottom=0.165, left = 0.16)
plt.savefig("/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_nuclear_data/violin_int_div_time_vs_lin_nuc_vol_growth_rate_r_%s_z_%s_%s.pdf"%(radius, zThres, region), format='pdf', dpi=300)
plt.close(fig)


cc_frac_new_mod = [cc_frac_new[i] for i in range(len(cc_frac_new)) if cc_nuc_volsU[i] != 0 and cc_nuc_volsD[i] != 0 ] 
cc_nucgr_new_mod = [cc_nuc_gr_new[i] for i in range(len(cc_frac_new)) if cc_nuc_volsU[i] != 0 and cc_nuc_volsD[i] != 0 ] 


CCBins = [0, 0.25, 0.5 , 0.75, 1.0]
grBins = [[] for k in range(len(CCBins)-1)]
meanNVGR = np.mean(cc_nucgr_new_mod)

for k in range(len(cc_frac_new_mod)):
    for l in range(len(CCBins)-1): 
        if (cc_frac_new_mod[k] >= CCBins[l] and cc_frac_new_mod[k] < CCBins[l + 1]): 
            grBins[l].append(cc_nucgr_new_mod[k])

labs = []
for t in range(len(CCBins)-1):
    labs.append(str(CCBins[t]) + '~' + str(CCBins[t+1]))

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
sm.graphics.violinplot(grBins, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs', 'violin_lw':0,
                                  'label_fontsize':'large',
                                  'label_rotation':30})
ax.set_xlabel(' frac. inter-division time ')
ax.set_ylabel(' nuc. vol. rel. growth rate ')
ax.set_ylim([-0.05,0.1])
fig.set_size_inches(5.0, 4.2)
plt.gcf().subplots_adjust(bottom=0.165, left = 0.16)
plt.savefig("/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_nuclear_data/violin_int_div_time_vs_rel_nuc_vol_growth_rate_r_%s_z_%s_%s.pdf"%(radius, zThres, region), format='pdf', dpi=300)
plt.close(fig)


CCBins = [0, 0.25, 0.5 , 0.75, 1.0]
grBins = [[] for k in range(len(CCBins)-1)]

for k in range(len(cc_frac_new)):
    for l in range(len(CCBins)-1): 
        if (cc_frac_new[k] >= CCBins[l] and cc_frac_new[k] < CCBins[l + 1]): 
            grBins[l].append(cc_vol_gr_new[k])

labs = []
for t in range(len(CCBins)-1):
    labs.append(str(CCBins[t]) + '~' + str(CCBins[t+1]))

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
sm.graphics.violinplot(grBins, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs', 'violin_lw':0,
                                  'label_fontsize':'large',
                                  'label_rotation':30})
ax.set_xlabel(' frac. inter-division time ')
ax.set_ylabel(' vol. rel. growth rate ')
ax.set_ylim([-0.05,0.1])
fig.set_size_inches(5.0, 4.2)
plt.gcf().subplots_adjust(bottom=0.165, left = 0.16)
plt.savefig("/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_nuclear_data/violin_int_div_time_vs_rel_vol_growth_rate_r_%s_z_%s_%s.pdf"%(radius, zThres, region), format='pdf', dpi=300)
plt.close(fig)


cc_nuc_gr = [(np.array(all_norm_nuc_vols[i][j][1:-1]) - np.array(all_norm_nuc_vols[i][j][:-2]))/(0.5*np.array(all_norm_nuc_vols[i][j][1:-1]) + 0.5*np.array(all_norm_nuc_vols[i][j][:-2]))*1.0/(np.array(all_times[i][j][1:-1])  - np.array(all_times[i][j][:-2])) for i in range(len(all_norm_nuc_vols)) for j in range(len(all_norm_nuc_vols[i])) ]
cc_nuc_volsU = [all_norm_nuc_vols[i][j][1:-1] for i in range(len(all_norm_nuc_vols)) for j in range(len(all_norm_nuc_vols[i]))]
cc_nuc_volsD = [all_norm_nuc_vols[i][j][:-2] for i in range(len(all_norm_nuc_vols)) for j in range(len(all_norm_nuc_vols[i]))]
	

cc_nuc_gr_new = [item for sublist in cc_nuc_gr for item in sublist] 
cc_nuc_volsU = [item for sublist in cc_nuc_volsU for item in sublist] 
cc_nuc_volsD = [item for sublist in cc_nuc_volsD for item in sublist]

print len(cc_frac_new)
print len(cc_nuc_gr_new)
cc_frac_new_mod = [cc_frac_new[i] for i in range(len(cc_frac_new)) if cc_nuc_volsU[i] != 0 and cc_nuc_volsD[i] != 0 ] 
cc_nucgr_new_mod = [cc_nuc_gr_new[i] for i in range(len(cc_frac_new)) if cc_nuc_volsU[i] != 0 and cc_nuc_volsD[i] != 0 ] 

CCBins = [0, 0.25, 0.5 , 0.75, 1.0]
grBins = [[] for k in range(len(CCBins)-1)]
meanNVGR = np.mean(cc_nucgr_new_mod)

for k in range(len(cc_frac_new_mod)):
    for l in range(len(CCBins)-1): 
        if (cc_frac_new_mod[k] >= CCBins[l] and cc_frac_new_mod[k] < CCBins[l + 1]): 
            grBins[l].append(cc_nucgr_new_mod[k])

labs = []
for t in range(len(CCBins)-1):
    labs.append(str(CCBins[t]) + '~' + str(CCBins[t+1]))

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
sm.graphics.violinplot(grBins, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs', 'violin_lw':0,
                                  'label_fontsize':'large',
                                  'label_rotation':30})
ax.set_xlabel(' frac. inter-division time')
ax.set_ylabel(' norm. nuc. vol. rel. growth rate')
ax.set_ylim([-0.05,0.1])
fig.set_size_inches(6.5, 5.5)
plt.gcf().subplots_adjust(bottom=0.165, left = 0.16)
plt.savefig("/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_nuclear_data/violin_int_div_time_vs_rel_norm_nuc_vol_growth_rate_r_%s_z_%s_%s.pdf"%(radius, zThres, region), format='pdf', dpi=300)
plt.close(fig)


CCBins = [0, 0.25, 0.5 , 0.75, 1.0]
nucVolBins = [[] for k in range(len(CCBins)-1)]

for k in range(len(cc_phase)):
    for l in range(len(CCBins)-1): 
        if (cc_phase[k] >= CCBins[l] and cc_phase[k] < CCBins[l + 1]): 
            nucVolBins[l].append(nuc_vol_cc_phase[k])
print nucVolBins
labs = []
for t in range(len(CCBins)-1):
    labs.append('[' + str(CCBins[t]) + ',' + str(CCBins[t+1]) + ']')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
sm.graphics.violinplot(nucVolBins, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs', 'violin_lw':0,
                                  'label_fontsize':'medium',
                                  'label_rotation':0})
ax.set_xlabel(' frac. inter-division time ')
ax.set_ylabel(' nuclear volume ')
ax.set_ylim([10, 95])
plt.yticks([20,40,60,80])
fig.set_size_inches(4.5, 2.5)
plt.gcf().subplots_adjust(bottom=0.25, left = 0.16)
plt.savefig("/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_nuclear_data/violin_int_div_time_vs_nuc_volume_r_%s_z_%s_%s.pdf"%(radius, zThres, region), format='pdf', dpi=300)
plt.close(fig)


nucVolBins = [[] for k in range(len(CCBins)-1)]

for k in range(len(cc_phase)):
    for l in range(len(CCBins)-1): 
        if (cc_phase[k] >= CCBins[l] and cc_phase[k] < CCBins[l + 1]): 
            nucVolBins[l].append(nuc_vol_cc_phase[k]/vol_cc_phase[k])
print nucVolBins
labs = []
for t in range(len(CCBins)-1):
    labs.append('['+str(CCBins[t]) + ',' + str(CCBins[t+1])+']')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
sm.graphics.violinplot(nucVolBins, labels = labs, ax=ax, 
                       plot_opts={'cutoff_val':5, 'cutoff_type':'abs', 'violin_lw':0,
                                  'label_fontsize':'medium',
                                  'label_rotation':0})
ax.set_xlabel(' frac. inter-division time ')
ax.set_ylabel(' nuclear volume/cell volume ')
ax.set_ylim([0.05, 0.55])
#plt.yticks([20,40,60,80])
fig.set_size_inches(4.5, 2.5)
plt.gcf().subplots_adjust(bottom=0.25, left = 0.16)
plt.savefig("/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/results_nuclear_data/violin_int_div_time_vs_frac_nuc_volume_r_%s_z_%s_%s.pdf"%(radius, zThres, region), format='pdf', dpi=300)
plt.close(fig)


