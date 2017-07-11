%%
% -------------------------------------------------------------------------
% Import lineage tracking data
% Input:
% pathToFiles_track - complete path to directory where all ALT parent label
% files are located
% Output:
% data_tracking - cell array with all the parent label tables (each element
% in the array corresponds to a time point
% -------------------------------------------------------------------------

function [data_tracking,filelist_track] = importDataTrackingALT(pathToFiles_track)

    filelist_track=dir(strcat(pathToFiles_track,'*.txt')); 

    for i=1:size(filelist_track,1)
        
        data_tracking{i} = makeALTtable(strcat(pathToFiles_track,filelist_track(i).name));
        
    end
    
end