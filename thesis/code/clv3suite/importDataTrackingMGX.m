%%
% -------------------------------------------------------------------------
% Import lineage tracking data
%
% Input:
% pathToFiles_track - complete path to directory where all MGX parent label
% files are located
%
% Output:
% data_tracking - cell array with all the parent label tables (each element
% in the array corresponds to a time point
% -------------------------------------------------------------------------

function [data_tracking] = importDataTrackingMGX(pathToFiles_track)

    filelist_track=dir(strcat(pathToFiles_track,'*.csv')); % Assumes Costanza data files are .txt

    for i=1:size(filelist_track,1)

        data_tracking{i} = csvread(strcat(pathToFiles_track,filelist_track(i).name),1);
        data_tracking{i}(:,[1,2])=data_tracking{i}(:,[2,1]);
    
    end
    
end