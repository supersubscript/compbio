#config file

rootPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant15/division_data/"
savePath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant15/statistics/"
workspacePath = "/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/workspace/"

#time, best ALT parameter, intensity image resolution, z-stretching factor computed from rapid scan,real time assuming dawn at time 0.
timesList = [
["0hrs", 3.0,  (0.2396150, 0.2396150, 0.26), 1.04, 10 + 56.0/60 -4], 
["4hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.105, 14 + 52.0/60 -4], 
["8hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.084, 18 + 46.0/60 -4], 
["12hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.105, 22 + 52.0/60 -4], 
["16hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.062, 24 + 2 + 39.0/60  - 4], 
["20hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.014, 24 + 6 + 33.0/60 -4],
["24hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.066, 24 + 10 + 49.0/60 -4], 
["28hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.066, 24 + 15 + 11.0/60 -4], 
["32hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.214, 24 + 19 + 26.0/60 -4], 
["36hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.300, 24 + 23 + 16.0/60 -4], 
["40hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.430, 48 + 3 + 23.0/60 -4], 
["44hrs", 5.0, (0.2196471, 0.2196471, 0.26), 1.582, 48 + 7 + 59.0/60 -4], 
["48hrs", 5.0, (0.2196471, 0.2196471, 0.26), 1.300, 48 + 10 + 27.0/60 -4], 
["52hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.539, 48 + 15 + 42.0/60 -4], 
["56hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.777, 48 + 19 + 23.0/60 -4], 
["60hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.517, 48 + 22 + 55.0/60 -4], 
["64hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.842, 72 + 3 + 9.0/60 -4],
["68hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.638, 72 + 7 + 17.0/60 -4], 
["72hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.742, 72 + 11 -4],
["76hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.474, 72 + 15 + 24.0/60 -4], 
["80hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.019, 72 + 18 + 39.0/60 -4],
["84hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.365, 72 + 22 + 53.0/60 -4]
]


altPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant15/tracking_data/"
ALTDataPathList =  [
altPath + 'matchingScoresListNew_0hrs_to_4hrs.pkl', 
altPath +'matchingScoresListNew_4hrs_to_8hrs.pkl', 
altPath +'matchingScoresListNew_8hrs_to_12hrs.pkl', 
altPath +'matchingScoresListNew_12hrs_to_16hrs.pkl', 
altPath +'matchingScoresListNew_16hrs_to_20hrs.pkl',
altPath +'matchingScoresListNew_20hrs_to_24hrs.pkl', 
altPath +'matchingScoresListNew_24hrs_to_28hrs.pkl', 
altPath +'matchingScoresListNew_28hrs_to_32hrs.pkl',
altPath +'matchingScoresListNew_32hrs_to_36hrs.pkl', 
altPath +'matchingScoresListNew_36hrs_to_40hrs.pkl', 
altPath +'matchingScoresListNew_40hrs_to_44hrs.pkl',
altPath +'matchingScoresListNew_44hrs_to_48hrs.pkl', 
altPath +'matchingScoresListNew_48hrs_to_52hrs.pkl',
altPath +'matchingScoresListNew_52hrs_to_56hrs.pkl',
altPath +'matchingScoresListNew_56hrs_to_60hrs.pkl',
altPath +'matchingScoresListNew_60hrs_to_64hrs.pkl', 
altPath +'matchingScoresListNew_64hrs_to_68hrs.pkl',
altPath +'matchingScoresListNew_68hrs_to_72hrs.pkl', 
altPath +'matchingScoresListNew_72hrs_to_76hrs.pkl', 
altPath +'matchingScoresListNew_76hrs_to_80hrs.pkl',
altPath +'matchingScoresListNew_80hrs_to_84hrs.pkl']

segPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant15/segmentation_data/"

segDataPathList =  [
segPath + '0hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '4hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '8hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '12hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '16hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '20hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '24hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '28hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '32hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '36hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '40hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '44hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '48hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '52hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '56hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '60hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '64hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '68hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '72hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '76hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '80hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl', 
segPath + '84hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3_top.pkl']



finalSegmentationPathList = [
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/0hrs/0hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/4hrs/4hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/8hrs/8hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/12hrs/12hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/16hrs/16hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/20hrs/20hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/24hrs/24hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/28hrs/28hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/32hrs/32hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/36hrs/36hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/40hrs/40hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/44hrs/44hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/48hrs/48hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/52hrs/52hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/56hrs/56hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/60hrs/60hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/64hrs/64hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/68hrs/68hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/72hrs/72hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/76hrs/76hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/80hrs/80hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/84hrs/84hrs_plant15_trim-acylYFP_hmin_2_asf_0_s_2.00_clean_3.tif"
]

segmentationParameters = [2, [0], 2.0] #[hmin, [asfList], sigma]
cleanParameter = "3" #radius

intensityImagePathList = [
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/0hrs/0hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/4hrs/4hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/8hrs/8hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/12hrs/12hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/16hrs/16hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/20hrs/20hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/24hrs/24hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/28hrs/28hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/32hrs/32hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/36hrs/36hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/40hrs/40hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/44hrs/44hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/48hrs/48hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/52hrs/52hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/56hrs/56hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/60hrs/60hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/64hrs/64hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/68hrs/68hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/72hrs/72hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/76hrs/76hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/80hrs/80hrs_plant15_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant15/plant15_topFiles_newAreas/84hrs/84hrs_plant15_trim-acylYFP.tif"
]

nucPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant15/nuclear_data/"
nuclearVolPathList = [
nucPath + "BOAvoxel_t00_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t04_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t08_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t08_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t16_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t20_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t24_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t28_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t32_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t36_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t40_volumes_segmentationLabels.pkl",
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
