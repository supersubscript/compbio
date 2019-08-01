

#paths and files
rootPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant4/division_data/"
workspacePath = "/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/workspace"
savePath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant4/statistics/"
cleanParameter = "3" #for cleaning
segmentationParameters = [2, [0,1], 1.5] #hmin, asf values, sigma


#time, best ALT parameter for corresponding interval, resolutions, stretching factor, real time assuming dawn at time 0.
timesList = [
 ["0hrs", 5.0, (0.2396150, 0.2396150, 0.26), 1.191, 10 + 8.0/60 - 4], 
 ["4hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.582,  14 + 16.0/60 - 4],
 ["8hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.928,  18+ 24.0/60 - 4], 
 ["12hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 2.21,  22 + 25.0/60 - 4],
 ["16hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.170, 24 + 02 + 11.0/60 - 4],
 ["20hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.274,  24 + 05 + 55.0/60 - 4], 
["24hrs", 5.0, (0.2196471, 0.2196471, 0.26), 1.222,  24 + 10 + 13.0/60 - 4],
["28hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.586,  24 + 14 + 16.0/60 - 4], 
["32hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.322,  24 + 18 + 22.0/60 - 4], 
["36hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.473,  24 + 22 + 30.0/60 - 4],
["40hrs", 5.0, (0.2196471, 0.2196471, 0.26), 1.170,  48 + 02 + 10.0/60 - 4],
["44hrs", 3.0, (0.2196471, 0.2196471, 0.29), 1.643,  48 + 06 + 53.0/60 - 4],
["48hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.196,  48 + 9 + 58.0/60 - 4],
["52hrs", 5.0, (0.2196471, 0.2196471, 0.26), 1.257,  48 + 14 + 34.0/60 - 4],
["56hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.192, 48 + 18 + 12.0/60 - 4], 
["60hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.278, 48 + 22 + 12.0/60 - 4],
["64hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.322,  72 +  02 + 32.0/60 - 4],
["68hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.352, 72 + 06 + 09.0/60 - 4], 
["72hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.170, 72 + 10 + 18.0/60 - 4], 
["76hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.127,  72 + 14 + 23.0/60 - 4]
]

altPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant4/tracking_data/"
ALTDataPathList =  [
altPath + 'matchingScoresListNew_0hrs_to_4hrs.pkl', 
altPath + 'matchingScoresListNew_4hrs_to_8hrs.pkl', 
altPath + 'matchingScoresListNew_8hrs_to_12hrs.pkl', 
altPath + 'matchingScoresListNew_12hrs_to_16hrs.pkl',
altPath +   'matchingScoresListNew_16hrs_to_20hrs.pkl',
altPath +  'matchingScoresListNew_20hrs_to_24hrs.pkl',
altPath +  'matchingScoresListNew_24hrs_to_28hrs.pkl',
altPath +  'matchingScoresListNew_28hrs_to_32hrs.pkl', 
altPath + 'matchingScoresListNew_32hrs_to_36hrs.pkl',
altPath +  'matchingScoresListNew_36hrs_to_40hrs.pkl',
altPath +  'matchingScoresListNew_40hrs_to_44hrs.pkl',
altPath +  'matchingScoresListNew_44hrs_to_48hrs.pkl',
altPath +  'matchingScoresListNew_48hrs_to_52hrs.pkl',
 altPath + 'matchingScoresListNew_52hrs_to_56hrs.pkl',
 altPath + 'matchingScoresListNew_56hrs_to_60hrs.pkl',
altPath +  'matchingScoresListNew_60hrs_to_64hrs.pkl',
altPath +  'matchingScoresListNew_64hrs_to_68hrs.pkl',
altPath +  'matchingScoresListNew_68hrs_to_72hrs.pkl',
altPath +  'matchingScoresListNew_72hrs_to_76hrs.pkl'
]



intensityImagePathList = [
 "/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/0hrs/0hrs_plant4_trim-acylYFP.tif",
 "/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/4hrs/4hrs_plant4_trim-acylYFP.tif",
 "/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/8hrs/8hrs_plant4_trim-acylYFP.tif",
 "/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/12hrs/12hrs_plant4_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/16hrs/16hrs_plant4_trim-acylYFP.tif",
 "/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/20hrs/20hrs_plant4_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/24hrs/24hrs_plant4_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/28hrs/28hrs_plant4_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/32hrs/32hrs_plant4_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/36hrs/36hrs_plant4_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/40hrs/40hrs_plant4_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/44hrs/44hrs_plant4_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/48hrs/48hrs_plant4_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/52hrs/52hrs_plant4_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/56hrs/56hrs_plant4_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/60hrs/60hrs_plant4_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/64hrs/64hrs_plant4_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/68hrs/68hrs_plant4_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/72hrs/72hrs_plant4_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/76hrs/76hrs_plant4_trim-acylYFP.tif"
]

finalSegmentationPathList = [
 "/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/0hrs/0hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3.tif", #potential for errors
 "/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/4hrs/4hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3.tif", #potential for errors
 "/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/8hrs/8hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3.tif", #potential for errors
 "/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/12hrs/12hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3.tif", #potential for errors
 "/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/16hrs/16hrs_plant4_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif", 
 "/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/20hrs/20hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/24hrs/24hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/28hrs/28hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/32hrs/32hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/36hrs/36hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/40hrs/40hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/44hrs/44hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/48hrs/48hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/52hrs/52hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3.tif",  
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/56hrs/56hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/60hrs/60hrs_plant4_trim-acylYFP_hmin_2_asf_1_s_2.50_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/64hrs/64hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/68hrs/68hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/72hrs/72hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant4/76hrs/76hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3.tif"
]

segPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant4/segmentation_data/"
segDataPathList = [segPath + '0hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3_top.pkl', 
 segPath + '4hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3_top.pkl', 
segPath + '8hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3_top.pkl', 
segPath + '12hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3_top.pkl', 
segPath + '16hrs_plant4_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl',
segPath + '20hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3_top.pkl', 
segPath + '24hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3_top.pkl', 
segPath + '28hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3_top.pkl',
segPath + '32hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3_top.pkl', 
segPath + '36hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3_top.pkl', 
segPath + '40hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3_top.pkl', 
segPath + '44hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3_top.pkl', 
segPath + '48hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3_top.pkl', 
segPath + '52hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3_top.pkl',
segPath + '56hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3_top.pkl', 
segPath + '60hrs_plant4_trim-acylYFP_hmin_2_asf_1_s_2.50_clean_3_top.pkl', 
segPath + '64hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3_top.pkl', 
segPath + '68hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3_top.pkl', 
segPath + '72hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3_top.pkl', 
segPath + '76hrs_plant4_trim-acylYFP_hmin_2_asf_0_s_2.50_clean_3_top.pkl'
]

nucPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant4/nuclear_data/"
nuclearVolPathList = [
nucPath + "BOAvoxel_t00_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t04_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t08_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t12_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t16_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t20_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t20_volumes_segmentationLabels.pkl", #error here
nucPath + "BOAvoxel_t28_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t32_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t36_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t40_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t44_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t48_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t52_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t56_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t60_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t64_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t68_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t72_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t76_volumes_segmentationLabels.pkl"
]

