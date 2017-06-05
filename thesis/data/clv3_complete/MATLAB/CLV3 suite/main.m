%% SECTION 1
% Import files
% Indicate full path of folder containing all files related to a given time
% course. Recommend one folder per data type: raw and tracking.
%
% NOTE: inside the folder, files should be named such that MATLAB imports them in the
% correct order automatically (e.g file_0h.tif may be imported after
% file_16h.tif, so it sould be named file_00h.tif instead)

% Nuclear segmentation files: segmentation tif files for each raw stack
pathToFiles_segmented='/media/Seagate Backup Plus Drive/CLV3 project/plant#2/segmented_MARS/'; %insert full path here
filelist_segmented=dir(strcat(pathToFiles_segmented,'*.tif')); % assumes .tif file extension

% Import tracking data: csv files with MorphoGraphX file format
% Requires importDataTrackingALT.m file
%[data_tracking] = importDataTrackingALT('/home/jose/home3/Projects/KF_nuclear_shape/150817-SUN1-timecourse-N18/tracking/');%insert full path here

% % Raw intensity files: single channel tif files of raw intensity signal
% pathToFiles_raw=strcat(path,'/source_data/raw/'); %insert full path here
% filelist_raw=dir(strcat(pathToFiles_raw,'*.tif')); % assumes .tif file extension

% Import Costanza data: txt files of Costanza quantification for each stack
% NOT: Costanza exports .xls files, change extension to txt
% Requires importDataCostanza.m file
[data_costanza] = importDataCostanza_realID('./tables_plant2/')%insert full path here

%%

% Get image resolutions for the time course

resolutions=[];

for k=1:size(filelist_segmented,1)
    
    fname=filelist_segmented(k).name;
    a=strcat(pathToFiles_segmented,fname);
    info= imfinfo(a); 
    resolutions(k,1)=info(1).Width;
    resolutions(k,2)=info(1).Height;
    info=[];
    
end

%%
nuclearData={};

for k=1%:size(filelist_segmented,1)
   
    k %counter to follow progress 
    
    fname=filelist_segmented(k).name;%segmented filename
    
    % Get number of slices in both raw and segmented stacks
    a=strcat(pathToFiles_segmented,fname);%full path
    info = imfinfo(a); 
    num_images = numel(info);%nr of slices in segmented stack
    
    % Get region properties for all slices of the segmented stack
    for i=1:num_images % loops through each segmented slice
        
        i
        
        img= imread(a,i);
        %for each slice get the properties of each segmented region
        %(corresponding to a nucleus)
        % Each element on 'allFrames' cell array contains the list of
        % properties for each segmented nucleus in that slice
        allFrames{1,i}=regionprops(img,'Area','BoundingBox','Centroid','ConvexHull','Eccentricity','MajorAxisLength','MinorAxisLength','Orientation','PixelList');
        nr_cells(1,i)=size(allFrames{1,i},1); %nr of nuclei segmented in each slice
        
        
    end
    
    total_cells(k)=max(nr_cells(1,:)); %total number of segmented nuclei over all slices
    
    % Loop through all slices for each nucleus, select largest slice and
    % calculate shape parameters
    for m=1%:size(allFrames{1,i},1) 
        
        nr=[1:1:total_cells(k)];

        for j=1:size(nr,2)
            
            count_j(j)=j;
            
            allAreas=[];
    
            for i=1:size(allFrames(m,:),2)
                
                count_i(j,i)=i;
                
                if size(allFrames{m,i},1)>=nr(j)

                    allAreas(i)=allFrames{m,i}(nr(j)).Area;
                       
                else
                    allAreas(i)=0;
        
                end

            end
            
            %allAreas_count(j,:)=allAreas;
            [maxVal,index]=max(allAreas);%maximum area annotation 
            
            if(max(allAreas)==0)
                
                % Shape parameter annotation
                nuclearData{k,j}.id=nr(j);                
                nuclearData{k,j}.PixelList=[0,0,0];
                nuclearData{k,j}.Area=0;
                nuclearData{k,j}.Eccentricity=0;
                nuclearData{k,j}.AspectRatio=0;
                
                
            else
                % Shape parameter annotation
                nuclearData{k,j}=allFrames{m,index}(nr(j));
                nuclearData{k,j}.id=nr(j);
                nuclearData{k,j}.AspectRatio=nuclearData{k,j}.BoundingBox(3) / nuclearData{k,j}.BoundingBox(4);
                nuclearData{k,j}.PixelList(:,3)=index;
                nuclearData{k,j}.AreaEccentricity=nuclearData{k,j}.Area*nuclearData{k,j}.Eccentricity;
                
            end
            
        end
        
        %allAreas=[];
    
    end
    
    %allFrames={};
    %nr_cells=[];
    
end

%%
% Get all the centroids and bounding boxes for each nucleus on each time point from Yassin's
% nuclear ID segmentations

centroids_mars={};
boundingBoxes={};

for i=1:size(nuclearData_plant2,1)
    
    for j=1:size(nuclearData_plant2,2)
        
        if isempty(nuclearData_plant2{i,j})==0
        
            centroids_mars{i}(j,:)=nuclearData_plant2{i,j}.Centroid;
            boundingBoxes{i}(j,:)=nuclearData_plant2{i,j}.BoundingBox;
            
        end
        
    end
    
end

%%
centroids_mars={};

for k=1:size(filelist_segmented,1)
   
    k %counter to follow progress 
    
    fname=filelist_segmented(k).name;%segmented filename
    
    % Get number of slices in both raw and segmented stacks
    a=strcat(pathToFiles_segmented,fname);%full path
    info = imfinfo(a); 
    num_images = numel(info);%nr of slices in segmented stack
    
    % Get region properties for all slices of the segmented stack
    for i=1:num_images % loops through each segmented slice
        
        i;
        
        img= imread(a,i);
        %for each slice get the properties of each segmented region
        %(corresponding to a nucleus)
        % Each element on 'allFrames' cell array contains the list of
        % properties for each segmented nucleus in that slice
        allFrames{1,i}=regionprops(img,'Area','BoundingBox','Centroid','ConvexHull','Eccentricity','MajorAxisLength','MinorAxisLength','Orientation','PixelList');
        nr_cells(1,i)=size(allFrames{1,i},1); %nr of nuclei segmented in each slice
        
        
    end
    
    total_cells(k)=max(nr_cells(1,:)); %total number of segmented nuclei over all slices
    
    % Loop through all slices for each nucleus, select largest slice and
    % calculate shape parameters
    for m=1%:size(allFrames{1,i},1) 
        
        nr=[1:1:total_cells(k)];

        for j=2:size(nr,2)
            
            %count_j(j)=j;
            
            allPixels=[NaN,NaN,NaN];
    
            for i=1:size(allFrames(m,:),2)
                i;
                %count_i(j,i)=i;
                
                if size(allFrames{m,i},1)>=nr(j)
                    
                    addPixels(:,1)=allFrames{m,i}(nr(j)).PixelList(:,1);
                    addPixels(:,2)=allFrames{m,i}(nr(j)).PixelList(:,2);
                    addPixels(:,3)=i;
                    
                    allPixels=cat(1,allPixels,addPixels);
                    
                end
                
                addPixels=[];
                       
            end
            
            centroids_mars{k}(j,1)=nanmean(allPixels(:,1));
            centroids_mars{k}(j,2)=nanmean(allPixels(:,2));
            centroids_mars{k}(j,3)=nanmean(allPixels(:,3));
            
            
        end
        

    
    end
    
    allFrames={};
    nr_cells=[];
    
end



%%

centroids_mars_corrected={};
centroids_costanza_corrected={};

for i=1:size(centroids_mars,2)
    
    for j=1:size(centroids_mars{i},1)
        
        centroids_mars_corrected{i}(:,1)=centroids_mars{i}(:,1).*pixelsizes(1,i);
        centroids_mars_corrected{i}(:,2)=centroids_mars{i}(:,2).*pixelsizes(2,i);
        centroids_mars_corrected{i}(:,3)=centroids_mars{i}(:,3).*pixelsizes(3,i);
        
    end
    
end

for i=1:size(data_costanza,2)
    
    for j=1:size(data_costanza{i},1)
        
        centroids_costanza_corrected{i}(:,1)=data_costanza{i}(:,2)./pixelsizes(1,i);
        centroids_costanza_corrected{i}(:,2)=data_costanza{i}(:,3)./pixelsizes(2,i);
        centroids_costanza_corrected{i}(:,3)=data_costanza{i}(:,4)./pixelsizes(3,i);
        
    end
    
end

for i=1:size(centroids_costanza_corrected,2)
    
    for j=1:size(centroids_costanza_corrected{i},1)

        xq=centroids_costanza_corrected{i}(j,1);
        yq=centroids_costanza_corrected{i}(j,2);
        zq=centroids_costanza_corrected{i}(j,3);
        
        for k=1:size(centroids_mars_corrected{i},1)

            xv=centroids_mars_corrected{i}(k,1);
            yv=centroids_mars_corrected{i}(k,2);
            zv=centroids_mars_corrected{i}(k,3);
            
            if round(xq)==round(xv) && round(yq)==round(yv) && round(zq)==round(zv)
                
                correspondences{i}(z,1)=k;
                correspondences{i}(z,2)=data_costanza{i}(j,1);
                z=z+1;
                
            end
            
        end
        
    end
    
    z=1;

end


%%
% Get centroids for mars and costanza nuclear segmentations, corrected for
% pixel sizes

centroids_mars_corrected={};
centroids_costanza_corrected={};

for i=1:size(centroids_mars,2)
    
    for j=1:size(centroids_mars{i},1)
        
        centroids_mars_corrected{i}(:,1)=centroids_mars{i}(:,1).*pixelsizes(1,i);
        centroids_mars_corrected{i}(:,2)=centroids_mars{i}(:,2).*pixelsizes(2,i);
        centroids_mars_corrected{i}(:,3)=centroids_mars{i}(:,3).*pixelsizes(3,i);
        
    end
    
end

for i=1:size(data_costanza,2)
    
    for j=1:size(data_costanza{i},1)
        
        centroids_costanza_corrected{i}(:,1)=data_costanza{i}(:,2)./pixelsizes(1,i);
        centroids_costanza_corrected{i}(:,2)=data_costanza{i}(:,3)./pixelsizes(2,i);
        centroids_costanza_corrected{i}(:,3)=data_costanza{i}(:,4)./pixelsizes(3,i);
        
    end
    
end

% Get correspondence tables for all nuclei in all time points by defining
% the closest pair of nuclei between costanza and mars segmentations

z=1;
for i=1:size(centroids_costanza_corrected,2)
    
    for j=1:size(centroids_costanza_corrected{i},1)
    
        xq=centroids_costanza_corrected{i}(j,1);
        yq=centroids_costanza_corrected{i}(j,2);
        zq=centroids_costanza_corrected{i}(j,3);
        
        distances=[];
        
        for k=1:size(centroids_mars_corrected{i},1)
            
            xv=centroids_mars_corrected{i}(k,1);
            yv=centroids_mars_corrected{i}(k,2);
            zv=centroids_mars_corrected{i}(k,3);
            
            X = [xq yq zq;xv yv zv];
            distances(k) = pdist(X,'euclidean');

        end
        
        [a,b]=min(distances);
        
        correspondences{i}(z,1)=b;
        correspondences{i}(z,2)=data_costanza{i}(j,1);
        z=z+1;
        
    end
    
    
    z=1;

end


%%
% Write correspondence tables
time=[0:4:80];
for i=1:size(correspondences,2)
    
    fname = sprintf('ID_correspondence_plant#2_%d.txt',time(i));
    csvwrite(fname,correspondences{i})
    
end

for i=1:size(centroids_mars_corrected,2)
    
    fname = sprintf('centroids_mars_plant#2_%d.txt',time(i));
    csvwrite(fname,centroids_mars_corrected{i})
    
end

for i=1:size(centroids_costanza_corrected,2)
    
    fname = sprintf('centroids_costanza_plant#2_%d.txt',time(i));
    csvwrite(fname,centroids_costanza_corrected{i})
    
end




%%
% % Check if costanza centroids are inside MARS bounding box
% 
% correspondences={};
% z=1;
% 
% for i=1:size(data_costanza,2)
%     
%     for j=1:size(data_costanza{i},1)
% 
%         xq=data_costanza{i}(j,2)/pixelsizes(1,i);
%         yq=data_costanza{i}(j,3)/pixelsizes(2,i);
%         
%         for k=2:size(boundingBoxes{i},1)
% 
%             xv=[boundingBoxes{i}(k,1)*pixelsizes(1,i),boundingBoxes{i}(k,1)*pixelsizes(1,i)+boundingBoxes{i}(k,3)*pixelsizes(1,i)];
%             yv=[boundingBoxes{i}(k,2)*pixelsizes(2,i),boundingBoxes{i}(k,2)*pixelsizes(2,i)+boundingBoxes{i}(k,4)*pixelsizes(2,i)];
% 
%             in = inpolygon(xq,yq,xv,yv);
%             
%             if in==1
%                 
%                 correspondences{i}(z,1)=k;
%                 correspondences{i}(z,2)=data_costanza{i}(j,1);
%                 z=z+1;
%                 
%             end
%             
%         end
%         
%     end
%     
%     z=1;
% 
% end
% 
% rectangle('Position', [nuclearData_plant2{1,2}.BoundingBox(1),nuclearData_plant2{1,2}.BoundingBox(2),nuclearData_plant2{1,2}.BoundingBox(3),nuclearData_plant2{1,2}.BoundingBox(4)],...
%   'EdgeColor','r','LineWidth',2 )
% 
% 
% 
