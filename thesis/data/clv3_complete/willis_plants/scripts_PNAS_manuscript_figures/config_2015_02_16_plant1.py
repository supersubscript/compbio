
#paths and files
rootPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant1/division_data/"
savePath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant1/statistics/"
radius = "3" #for cleaning


#time, best ALT parameter for corresponding interval, resolutions, stretching factor, real time assuming dawn at time 0.


timesList = [
["0hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.06, 10 + 5.0/60 - 4],
#["4hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 2.77, 13.54 - 4],
["8hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.47, 18.0 - 4],
["12hrs", 3.0,  (0.2196471, 0.2196471, 0.26), 1.23, 21 + 57.0/60 - 4],
["16hrs", 5.0,  (0.2196471, 0.2196471, 0.26), 1.17, 24 + 1 + 37.0/60 - 4],
["20hrs", 5.0, (0.2396150, 0.2396150, 0.26), 1.20, 24 + 5 + 31.0/60 - 4],
["24hrs", 3.0, (0.2396150, 0.2396150, 0.26), 1.33, 24 + 9 + 57.0/60 - 4],
["28hrs", 3.0, (0.2396150, 0.2396150, 0.26), 1.33, 24 + 13 + 43.0/60 - 4], 
["32hrs", 3.0, (0.2396150, 0.2396150, 0.26), 1.37, 24 + 17 + 51.0/60 - 4],
["36hrs", 5.0, (0.2396150, 0.2396150, 0.26), 1.26, 24 + 21 + 48.0/60 - 4],
["40hrs", 5.0, (0.2396150, 0.2396150, 0.26), 1.23, 48 + 1 + 46.0/60 - 4],
["44hrs", 5.0, (0.2396150, 0.2396150, 0.26), 1.23, 48 + 5 + 56.0/60 - 4],
["48hrs", 3.0, (0.2396150, 0.2396150, 0.26), 1.22, 48 + 9 + 47.0/60 - 4],
["52hrs", 3.0, (0.2396150, 0.2396150, 0.26), 1.00, 48 + 13 + 56.0/60 - 4],
["56hrs", 3.0, (0.2396150, 0.2396150, 0.26), 1.15, 48 + 18 + 08.0/60 - 4],
["60hrs", 3.0, (0.2396150, 0.2396150, 0.26), 1.28, 48 + 21 + 48.0/60 - 4],
["64hrs", 3.0, (0.2396150, 0.2396150, 0.26), 1.36, 72 + 1 + 49.0/60 - 4],
["68hrs", 2.5, (0.2396150, 0.2396150, 0.26), 1.30, 72 + 5 + 44.0/60 - 4],
["72hrs", 3.0, (0.2396150, 0.2396150, 0.26), 1.12, 72 + 10 + 24.0/60 - 4],
["76hrs", 3.0, (0.2396150, 0.2396150, 0.26), 1.06, 72 + 13 + 57.0/60 - 4]]

altPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant1/tracking_data/"
ALTDataPathList = [
altPath + "matchingScoresListNew_0hrs_to_8hrs.pkl",
#altPath + "matchingScoresListNew_0hrs_to_4hrs.pkl",
#altPath + "matchingScoresListNew_4hrs_to_8hrs.pkl",
altPath + "matchingScoresListNew_8hrs_to_12hrs.pkl",
altPath + "matchingScoresListNew_12hrs_to_16hrs.pkl",
altPath + "matchingScoresListNew_16hrs_to_20hrs.pkl",
altPath + "matchingScoresListNew_20hrs_to_24hrs.pkl",
altPath + "matchingScoresListNew_24hrs_to_28hrs.pkl",
altPath + "matchingScoresListNew_28hrs_to_32hrs.pkl",
altPath + "matchingScoresListNew_32hrs_to_36hrs.pkl",
altPath + "matchingScoresListNew_36hrs_to_40hrs.pkl",
altPath + "matchingScoresListNew_40hrs_to_44hrs.pkl",
altPath + "matchingScoresListNew_44hrs_to_48hrs.pkl",
#altPath + "matchingScoresListNew_48hrs_to_56hrs.pkl",
altPath + "matchingScoresListNew_48hrs_to_52hrs.pkl",
altPath + "matchingScoresListNew_52hrs_to_56hrs.pkl",
altPath + "matchingScoresListNew_56hrs_to_60hrs.pkl",
altPath + "matchingScoresListNew_60hrs_to_64hrs.pkl",
altPath + "matchingScoresListNew_64hrs_to_68hrs.pkl",
altPath + "matchingScoresListNew_68hrs_to_72hrs.pkl",
altPath + "matchingScoresListNew_72hrs_to_76hrs.pkl",
]

segPath = "/home/lisa/home3/2015_02_16_NPA_PGMYCR_timecourse/cell_size_growth_regulation_all_final_data/plant1/segmentation_data/"
segDataPathList = [
segPath + "0hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl",
#segPath + "4hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl",
segPath + "8hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl",
segPath + "12hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl",
segPath + "16hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl",
segPath + "20hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl",
segPath + "24hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl",
segPath + "28hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl",
segPath + "32hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl",
segPath + "36hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl",
segPath + "40hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl",
segPath + "44hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl",
segPath + "48hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl",
segPath + "52hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl",
segPath + "56hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl",
segPath + "60hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl",
segPath + "64hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl",
segPath + "68hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl",
segPath + "72hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl",
segPath + "76hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3_top.pkl"
]


finalSegmentationPathList = [
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/0hrs/0hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/4hrs/4hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/8hrs/8hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/12hrs/12hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/16hrs/16hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/20hrs/20hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/24hrs/24hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/28hrs/28hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/32hrs/32hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/36hrs/36hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/40hrs/40hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/44hrs/44hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/48hrs/48hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/52hrs/52hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/56hrs/56hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/60hrs/60hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/64hrs/64hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/68hrs/68hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/72hrs/72hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/76hrs/76hrs_plant1_trim-acylYFP_hmin_2_asf_1_s_2.00_clean_3.tif"
]

intensityImagePathList = [
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/0hrs/0hrs_plant1_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/4hrs/4hrs_plant1_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/8hrs/8hrs_plant1_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/12hrs/12hrs_plant1_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/16hrs/16hrs_plant1_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/20hrs/20hrs_plant1_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/24hrs/24hrs_plant1_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/28hrs/28hrs_plant1_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/32hrs/32hrs_plant1_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/36hrs/36hrs_plant1_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/40hrs/40hrs_plant1_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/44hrs/44hrs_plant1_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/48hrs/48hrs_plant1_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/52hrs/52hrs_plant1_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/56hrs/56hrs_plant1_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/60hrs/60hrs_plant1_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/64hrs/64hrs_plant1_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/68hrs/68hrs_plant1_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/72hrs/72hrs_plant1_trim-acylYFP.tif",
"/Users/lisawillis/NPA_data_analysis/seg_data/plant1/76hrs/76hrs_plant1_trim-acylYFP.tif"
]

