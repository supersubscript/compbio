
#paths and files
rootPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant18/division_data/"
workspacePath = "/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/workspace"
savePath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant18/statistics/"
cleanParameter = "3" #for cleaning
segmentationParameters = [2, [0],3.0] #hmin, asf values, sigma

#time, best ALT parameter for corresponding interval, resolutions, stretching factor, real time assuming dawn at time 0.
timesList = [
["0hrs", 3.0, (0.2396150, 0.2396150, 0.26), 0.997, 11 + 28.0/60 - 4], 
["4hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 0.975,  15 + 15.0/60 - 4],
["8hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.062,  18 + 59.0/60 - 4],
["12hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.018,  23 + 18.0/60 - 4],
["16hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.062, 24 + 02 + 53.0/60 - 4],
["20hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.014,  24 + 06 + 43.0/60 - 4], 
["24hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.066,  24 + 11 + 00.0/60 - 4],
["28hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.014,  24 + 15 + 04.0/60 - 4], 
["32hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.235,  24 + 20 + 02.0/60 - 4],
["36hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.257,  24 + 23 + 27.0/60 - 4], 
#["40hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.278,  48 + 03 + 19.0/60 - 4], 
["44hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.322,  48 + 8 + 29.0/60 - 4],
["48hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.404,  48 + 10 + 38.0/60 - 4],
["52hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.538,  48 + 16 + 20.0/60 - 4],
["56hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.582, 48 + 19 + 23.0/60 - 4],
["60hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.408, 48 + 23 + 06.0/60 - 4],
["64hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.300,  72 + 03 + 22.0/60 - 4],
["68hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.56, 72 + 07 + 05.0/60 - 4],
["72hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.222, 72 + 11 + 10.0/60 - 4],
["76hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.235,  72 + 15 + 41.0/60 - 4],
["80hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.040, 72 + 18 + 53.0/60 - 4],
["84hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.495,  72 + 23 + 07.0/60 - 4]
]

intensityImagePathList = [
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/0hrs/0hrs_plant18_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/4hrs/4hrs_plant18_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/8hrs/8hrs_plant18_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/12hrs/12hrs_plant18_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/16hrs/16hrs_plant18_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/20hrs/20hrs_plant18_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/24hrs/24hrs_plant18_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/28hrs/28hrs_plant18_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/32hrs/32hrs_plant18_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/36hrs/36hrs_plant18_trim-acylYFP.tif", #mapping 36--> 44 hrs at 3.0 works
#"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/40hrs/40hrs_plant18_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/44hrs/44hrs_plant18_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/48hrs/48hrs_plant18_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/52hrs/52hrs_plant18_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/56hrs/56hrs_plant18_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/60hrs/60hrs_plant18_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/64hrs/64hrs_plant18_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/68hrs/68hrs_plant18_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/72hrs/72hrs_plant18_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/6hrs/76hrs_plant18_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/80hrs/80hrs_plant18_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/84hrs/84hrs_plant18_trim-acylYFP.tif"
]


finalSegmentationPathList = [
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/0hrs/0hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/4hrs/4hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/8hrs/8hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif", #missing L1 cell
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/12hrs/12hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/16hrs/16hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/20hrs/20hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/24hrs/24hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/28hrs/28hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/32hrs/32hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/36hrs/36hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif", 
#"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/40hrs/40hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/44hrs/44hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/48hrs/48hrs_plant18_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/52hrs/52hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",  
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/56hrs/56hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/60hrs/60hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/64hrs/64hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/68hrs/68hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.50_clean_3.tif", #s = 2.50 possibly introduces changes in volume 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/72hrs/72hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/76hrs/76hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/80hrs/80hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant18/plant18_topFiles_newAreas/84hrs/84hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif"
]

altPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant18/tracking_data/"
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
altPath + 'matchingScoresListNew_36hrs_to_44hrs.pkl', 
altPath + 'matchingScoresListNew_44hrs_to_48hrs.pkl', 
altPath + 'matchingScoresListNew_48hrs_to_52hrs.pkl', 
altPath + 'matchingScoresListNew_52hrs_to_56hrs.pkl', 
altPath + 'matchingScoresListNew_56hrs_to_60hrs.pkl', 
altPath + 'matchingScoresListNew_60hrs_to_64hrs.pkl', 
altPath + 'matchingScoresListNew_64hrs_to_68hrs.pkl', 
altPath + 'matchingScoresListNew_68hrs_to_72hrs.pkl', 
altPath + 'matchingScoresListNew_72hrs_to_76hrs.pkl', 
altPath + 'matchingScoresListNew_76hrs_to_80hrs.pkl', 
altPath + 'matchingScoresListNew_80hrs_to_84hrs.pkl'
]

segPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant18/segmentation_data/"
segDataPathList = [
segPath + '0hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl', 
segPath +'4hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl', 
segPath +'8hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl', 
segPath +'12hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl', 
segPath +'16hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl', 
segPath +'20hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl', 
segPath +'24hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl', 
segPath +'28hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl', 
segPath +'32hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl', 
segPath +'36hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl', 
segPath +'44hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl', 
segPath +'48hrs_plant18_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath +'52hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl', 
segPath +'56hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl', 
segPath +'60hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl', 
segPath +'64hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl',
segPath +'68hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.50_clean_3_top.pkl', 
segPath +'72hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl', 
segPath +'76hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl', 
segPath +'80hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl',
segPath +'84hrs_plant18_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl'
]

nucPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant18/nuclear_data/"
nuclearVolPathList = [nucPath + "BOAvoxel_t00_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t04_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t08_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t12_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t16_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t20_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t24_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t28_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t32_volumes_segmentationLabels.pkl",
nucPath + "BOAvoxel_t36_volumes_segmentationLabels.pkl",
#nucPath + "BOAvoxel_t40_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t44_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t48_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t52_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t56_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t60_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t64_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t68_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t72_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t76_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t80_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t84_volumes_segmentationLabels.pkl"
]

