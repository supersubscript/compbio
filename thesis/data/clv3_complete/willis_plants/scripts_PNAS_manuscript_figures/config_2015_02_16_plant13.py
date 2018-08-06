
#paths and files
rootPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant13/division_data/"
workspacePath = "/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/workspace"

savePath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant13/statistics/"
cleanParameter = "3" #for cleaning
segmentationParameters = [2, [1], 2.5] #hmin, asf values, sigma

#time, best ALT parameter for corresponding interval, resolutions, stretching factor, real time assuming dawn at time 0, max radius in microns to include all cells
timesList = [
["0hrs", 3.0, (0.2396150, 0.2396150, 0.26), 1.062, 10 + 45.0/60 - 4, 40],
["4hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.322,  14 + 44.0/60 - 4],
["8hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.083,  18+ 37.0/60 - 4],
["12hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.083,  23 + 02.0/60 - 4],
["16hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.040, 24 + 02+ 26.0/60 - 4],
["20hrs", 5.0, (0.2196471, 0.2196471, 0.26), 1.040,  24 + 06 + 23.0/60 - 4], 
["24hrs", 5.0, (0.2196471, 0.2196471, 0.26), 1.040,  24 + 10 + 58.0/60 - 4],
["28hrs", 5.0, (0.2196471, 0.2196471, 0.26), 1.066,  24 + 14 + 43.0/60 - 4], 
["32hrs", 5.0, (0.2196471, 0.2196471, 0.26), 1.192,  24 + 19 + 15.0/60 - 4],
["36hrs", 5.0, (0.2196471, 0.2196471, 0.26), 1.538,  24 + 23 + 04.0/60 - 4],
["40hrs", 5.0, (0.2196471, 0.2196471, 0.26), 1.321,  48 + 03 + 13.0/60 - 4],
["44hrs", 5.0, (0.2196471, 0.2196471, 0.26), 1.235,  48 + 07 + 42.0/60 - 4],
["48hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.274,  48 + 10 + 13.0/60 - 4],
["52hrs", 5.0, (0.2196471, 0.2196471, 0.26), 1.43,  48 + 15 + 18.0/60 - 4],
["56hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.798, 48 + 18 + 56.0/60 - 4],
["60hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.798, 48 + 22 + 46.0/60 - 4],
["64hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.538,  72 +  02 + 59.0/60 - 4],
["68hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.456, 72 + 06 + 31.0/60 - 4],
["72hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.872, 72 + 10 + 46.0/60 - 4],
["76hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.387,  72 + 15 + 5.0/60 - 4],
["80hrs", 3.0, (0.2196471, 0.2196471, 0.26), 1.062,  72 + 19 + 18.0/60 - 4, 45],
["84hrs", 3.0, (0.2196471, 0.2196471, 0.26), 2.318,  72 + 22 + 44.0/60 - 4]
]

intensityImagePathList = [
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/0hrs/0hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/4hrs/4hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/8hrs/8hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/12hrs/12hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/16hrs/16hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/20hrs/20hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/24hrs/24hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/28hrs/28hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/32hrs/32hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/36hrs/36hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/40hrs/40hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/44hrs/44hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/48hrs/48hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/52hrs/52hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/56hrs/56hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/60hrs/60hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/64hrs/64hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/68hrs/68hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/72hrs/72hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/76hrs/76hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/80hrs/80hrs_plant13_trim-acylYFP.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/84hrs/84hrs_plant13_trim-acylYFP.tif"
]

finalSegmentationPathList = [
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/0hrs/0hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/4hrs/4hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/8hrs/8hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/12hrs/12hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/16hrs/16hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/20hrs/20hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/24hrs/24hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/28hrs/28hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/32hrs/32hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_2.50_clean_3.tif", #1.50 MAY BE BETTER BUT INTRODUCES l1 ERROR 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/36hrs/36hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/40hrs/40hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_2.50_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/44hrs/44hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/48hrs/48hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/52hrs/52hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_2.50_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/56hrs/56hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/60hrs/60hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_2.50_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/64hrs/64hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_2.50_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/68hrs/68hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif",
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/72hrs/72hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_2.50_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/76hrs/76hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/80hrs/80hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3.tif", 
"/home/lisa/Desktop/NPA_data_analysis/seg_data/plant13/84hrs/84hrs_plant13_trim-acylYFP_improved_hmin_2_asf_1_s_2.50_clean_3.tif"
]

altPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant13/tracking_data/"
ALTDataPathList =  [
altPath + 'matchingScoresListNew_0hrs_to_4hrs.pkl', 
altPath +'matchingScoresListNew_4hrs_to_8hrs.pkl', 
altPath +'matchingScoresListNew_8hrs_to_12hrs.pkl', 
altPath +'matchingScoresListNew_12hrs_to_16hrs.pkl',
altPath + 'matchingScoresListNew_16hrs_to_20hrs.pkl',
 altPath +'matchingScoresListNew_20hrs_to_24hrs.pkl',
altPath + 'matchingScoresListNew_24hrs_to_28hrs.pkl',
altPath + 'matchingScoresListNew_28hrs_to_32hrs.pkl', 
altPath +'matchingScoresListNew_32hrs_to_36hrs.pkl',
altPath + 'matchingScoresListNew_36hrs_to_40hrs.pkl',
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
 altPath +'matchingScoresListNew_80hrs_to_84hrs.pkl'
]

segPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant13/segmentation_data/"
segDataPathList =  [
segPath + '0hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl', 
segPath +'4hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl', 
segPath +'8hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl', 
segPath +'12hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl',
segPath + '16hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl', 
segPath +'20hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl', 
segPath +'24hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl',
segPath + '28hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl',
segPath + '32hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_2.50_clean_3_top.pkl', 
segPath +'36hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl', 
segPath +'40hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_2.50_clean_3_top.pkl', 
segPath +'44hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl', 
segPath +'48hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl', 
segPath +'52hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_2.50_clean_3_top.pkl', 
segPath +'56hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl', 
segPath +'60hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_2.50_clean_3_top.pkl', 
segPath +'64hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_2.50_clean_3_top.pkl', 
segPath +'68hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl', 
segPath +'72hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_2.50_clean_3_top.pkl', 
segPath +'76hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl', 
segPath +'80hrs_plant13_trim-acylYFP_hmin_2_asf_1_s_1.50_clean_3_top.pkl',
segPath +"84hrs_plant13_trim-acylYFP_improved_hmin_2_asf_1_s_2.50_clean_3_top.pkl"
]

nucPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant13/nuclear_data/"
nuclearVolPathList = [
nucPath + "BOAvoxel_t00_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t04_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t08_volumes_segmentationLabels.pkl",
nucPath +"BOAvoxel_t12_volumes_segmentationLabels.pkl",
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
