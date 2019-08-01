
#paths and files
rootPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant2/division_data/"
workspacePath = "/home/lisa/Desktop/NPA_data_analysis/seg_data/plant2/workspace"
savePath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant2/statistics/"
cleanParameter = "3" #for cleaning
segmentationParameters = [2, [1], 2.2] #hmin, asf values, sigma

#time, best ALT parameter for corresponding interval, resolutions, stretching factor, real time assuming dawn at time 0.
timesList = [
["0hrs", 3.0, (0.2396150, 0.2396150, 0.26), 1.408, 9 + 41.0/60 - 4],
["4hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.17,  14 + 02.0/60 - 4],
["8hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.127,  18 + 11.0/60 - 4],  
["12hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.105,  21 + 50.0/60 - 4],
["16hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.083, 24 + 01 + 48.0/60 - 4],
["20hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.092,  24 + 06 + 04.0/60 - 4], 
["24hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.43,  24 + 10 + 15.0/60 - 4],
["28hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.404,  24 + 13 + 57.0/60 - 4], 
["32hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.32,  24 + 18 + 05.0/60 - 4],
["36hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.19,  24 + 22 + 00.0/60 - 4],
["40hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.32,  48 + 02 + 19.0/60 - 4],
["44hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.13,  48 + 06 + 40.0/60 - 4], #repeat track --manual
["48hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.222,  48 + 11 + 02.0/60 - 4],
["52hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.127,  48 + 14 + 32.0/60 - 4], #repeat track --manual # 1 possible error in seg
["56hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.278, 48 + 18 + 00.0/60 - 4],
["60hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.3, 48 + 22 + 15.0/60 - 4],
["64hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.365,  72 + 02 + 05.0/60 - 4],
["68hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.274, 72 + 05 + 54.0/60 - 4],
["72hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.404, 72 + 9 + 52.0/60 - 4],
["76hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.083,  72 + 14 + 11.0/60 - 4]
]


intensityImagePathList = [
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/0hrs/0hrs_plant2_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/4hrs/4hrs_plant2_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/8hrs/8hrs_plant2_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/12hrs/12hrs_plant2_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/16hrs/16hrs_plant2_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/20hrs/20hrs_plant2_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/24hrs/24hrs_plant2_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/28hrs/28hrs_plant2_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/32hrs/32hrs_plant2_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/36hrs/36hrs_plant2_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/40hrs/40hrs_plant2_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/44hrs/44hrs_plant2_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/48hrs/48hrs_plant2_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/52hrs/52hrs_plant2_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/56hrs/56hrs_plant2_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/60hrs/60hrs_plant2_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/64hrs/64hrs_plant2_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/68hrs/68hrs_plant2_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/72hrs/72hrs_plant2_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/76hrs/76hrs_plant2_trim-acylYFP.tif"
]

altPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant2/tracking_data/"

ALTDataPathList =  [
altPath + 'matchingScoresListNew_0hrs_to_4hrs.pkl', 
altPath + 'matchingScoresListNew_4hrs_to_8hrs.pkl', 
altPath + 'matchingScoresListNew_8hrs_to_12hrs.pkl', 
altPath + 'matchingScoresListNew_12hrs_to_16hrs.pkl', 
altPath + 'matchingScoresListNew_16hrs_to_20hrs.pkl', 
altPath + 'matchingScoresListNew_20hrs_to_24hrs.pkl', 
altPath + 'matchingScoresListNew_24hrs_to_28hrs.pkl', 
altPath + 'matchingScoresListNew_28hrs_to_32hrs.pkl', 
altPath + 'matchingScoresListNew_32hrs_to_36hrs.pkl', 
altPath + 'matchingScoresListNew_36hrs_to_40hrs.pkl', 
altPath + 'matchingScoresListNew_40hrs_to_44hrs.pkl', 
altPath + 'matchingScoresListNew_44hrs_to_48hrs.pkl', 
altPath + 'matchingScoresListNew_48hrs_to_52hrs.pkl', 
altPath + 'matchingScoresListNew_52hrs_to_56hrs.pkl', 
altPath + 'matchingScoresListNew_56hrs_to_60hrs.pkl', 
altPath + 'matchingScoresListNew_60hrs_to_64hrs.pkl', 
altPath + 'matchingScoresListNew_64hrs_to_68hrs.pkl', 
altPath + 'matchingScoresListNew_68hrs_to_72hrs.pkl', 
altPath + 'matchingScoresListNew_72hrs_to_76hrs.pkl'
]


finalSegmentationPathList = [
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/0hrs/0hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif", 
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/4hrs/4hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif", 
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/8hrs/8hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_2.20_clean_3.tif", #1.5 gives 1 error
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/12hrs/12hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif", 
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/16hrs/16hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif", 
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/20hrs/20hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif", 
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/24hrs/24hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif", 
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/28hrs/28hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_2.20_clean_3.tif", #1.5 gives 2 errors 
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/32hrs/32hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif", 
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/36hrs/36hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif", 
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/40hrs/40hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif", 
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/44hrs/44hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif", 
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/48hrs/48hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif", 
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/52hrs/52hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif",  
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/56hrs/56hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif", 
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/60hrs/60hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/64hrs/64hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/68hrs/68hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/72hrs/72hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant2/76hrs/76hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif"
]

segPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant2/segmentation_data/"

segDataPathList =  [
segPath + "0hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl", 
segPath + "4hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl", 
segPath + "8hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_2.20_clean_3_top.pkl", #1.5 gives 1 error
segPath + "12hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl", 
segPath + "16hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl", 
segPath + "20hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl", 
segPath + "24hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl", 
segPath + "28hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_2.20_clean_3_top.pkl", #1.5 gives 2 errors 
segPath + "32hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl", 
segPath + "36hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl", 
segPath + "40hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl", 
segPath + "44hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl", 
segPath + "48hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl", 
segPath + "52hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl",  
segPath + "56hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl", 
segPath + "60hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl",
segPath + "64hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl",
segPath + "68hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl",
segPath + "72hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl",
segPath + "76hrs_plant2_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl"
]

nucPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant2/nuclear_data/"
nuclearVolPathList = [
nucPath + "t0_BOA_volumes_segmentationLabels.pkl",
nucPath + "t4_BOA_volumes_segmentationLabels.pkl",
nucPath + "t8_BOA_volumes_segmentationLabels.pkl",
nucPath + "t12_BOA_volumes_segmentationLabels.pkl",
nucPath + "t16_BOA_volumes_segmentationLabels.pkl",
nucPath + "t20_BOA_volumes_segmentationLabels.pkl",
nucPath + "t24_BOA_volumes_segmentationLabels.pkl",
nucPath + "t28_BOA_volumes_segmentationLabels.pkl",
nucPath + "t32_BOA_volumes_segmentationLabels.pkl",
nucPath + "t36_BOA_volumes_segmentationLabels.pkl",
nucPath + "t40_BOA_volumes_segmentationLabels.pkl",
nucPath + "t44_BOA_volumes_segmentationLabels.pkl",
nucPath + "t48_BOA_volumes_segmentationLabels.pkl",
nucPath + "t52_BOA_volumes_segmentationLabels.pkl",
nucPath + "t56_BOA_volumes_segmentationLabels.pkl",
nucPath + "t60_BOA_volumes_segmentationLabels.pkl",
nucPath + "t64_BOA_volumes_segmentationLabels.pkl",
nucPath + "t68_BOA_volumes_segmentationLabels.pkl",
nucPath + "t72_BOA_volumes_segmentationLabels.pkl",
nucPath + "t76_BOA_volumes_segmentationLabels.pkl"
]