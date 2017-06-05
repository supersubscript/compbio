Article: Willis, L., Refahi, Y., Wightman, R., Landrein, B., Teles, J., Huang, K. C., Meyerowitz E. M.  & Jönsson, H. (2016). Cell size and growth regulation in the Arabidopsis thaliana apical stem cell niche. Proceedings of the National Academy of Sciences, 113(51), E8238-E8246.


Data Folders:
>>plantX (Relabelled within article: plant 1--> plant 1, plant 2--> plant 2, plant 4 --> plant 3, plant 13 --> plant 4, plant 15 --> plant 5, plant 18 --> plant 6.)


>>plantX/processed_tiffs
Timecourses of membrane data (acylYFP) and CLAVATA3 data (clv3) after minimal processing, as described in "Supplementary Methods and Materials". This is prior to re-computation of z-resolution according to Table S10.  

>>plantX/segmentation_tiffs
Segmented processed tiffs using membrane data. Each cell is identified by an unique integer label. The background is labelled "1".

>>plantX/nuclear_segmentation_tiffs
Nuclear segmentations using processed tiffs for CLAVATA3 data. Each cell is identified by an unique integer label. The background is labelled "1".
The nuclei have the same label as the cells (in segmentation_tiffs files) in which they are. 
Plants 1, and 18 do not have Nuclear segmentations.

>>plantX/plantX_tiff_resolutions.txt
The correct (x,y,z)-resolutions for each tiff associated with each point in the timecourse, and the computed z-stretching factor as reported in Table S10.
To calculate cell properties (e.g. volume, surface area), divide z-resolution by the z-stretching factor.

>>plantX/tracking_data
Cell lineage data. timePointT_to_timePointT+4.txt : text file whose lines are mother cell labels in segmented image at time T : daughter cells in segmented image at time T + 4





