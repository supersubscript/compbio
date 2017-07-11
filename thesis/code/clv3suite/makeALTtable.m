% -------------------------------------------------------------------------
% Make lineage tracking data starting from ALT format
% To be used by function importDataTrackingALT.m
% -------------------------------------------------------------------------

function [mgx_table]=makeALTtable(filename, outFile)

    mgx_table=[];
    k=1;

    fid = fopen(filename,'rt');

    tline = fgetl(fid);

    while tline~=-1
   
        tline(strfind(tline, ':')) = [];
        tline(strfind(tline, '[')) = [];
        tline(strfind(tline, ']')) = [];
        tline(strfind(tline, ',')) = [];
    
        temp_vector=sscanf(tline,'%f');
    
        for i=2:size(temp_vector,1)
    
            mgx_table(k,1)=temp_vector(1);
            mgx_table(k,2)=temp_vector(i);
            k=k+1;
    
        end
    
        % temp_vector=[];
        tline = fgetl(fid);
    end

    fclose(fid);
    
    %csvwrite(outFile,mgx_table)

end