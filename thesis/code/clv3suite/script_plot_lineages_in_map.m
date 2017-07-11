% ----------------------------------------------
% Plot lineages in 2D map
% ----------------------------------------------

%%
% Import tracking from ALT formated files 
[data_tracking] = importDataTrackingALT('../../data/clv3_complete/PNAS/plant2/tracking_data/')


%%
L1_nuclei=[28:1:68];
nuclei_subset=[60 61];
time=[0:4:76]; %time vector (change sampling and end point)

% ------------------------------------
% Variables for colormap definition
% ------------------------------------

%variable=allConcentration_final;% select variable
%gradation=1000;%Defines color range
%maxVariable=max(variable(:));%Maximum value for variable
%minVariable=min(variable(:));%Minimum value for variable
%gradationLevel = round(maxVariable*gradation)-(round(minVariable*gradation));%Define color gradient
%colorVector = jet(gradationLevel);%Choose colormap
%colormap(jet)%Set colormap

for i=2:size(nuclearData,1) %Loop through each time point
    i
    % Get the L1 nuclei from the ALT tracking tables and check
    % correspondence for Costanza

    h = figure(i)
    hold

        for j=2:size(nuclearData,2)%Loop through each nucleus in a given time point
            j %counter to follow progress (nucleus number)
            %if isempty(nuclearData{i,j})==0 %if element is not empty
            if any(j==L1_nuclei)==1 %if element is not empty
                % Get centroid and orientation and calculate X and Y major
                % and minor axes
                centroid = nuclearData{i,j}.Centroid;
                orientation = nuclearData{i,j}.Orientation;
                xMajor=centroid(1) + [-1 1]*(nuclearData{i,j}.MajorAxisLength/2)*cosd(orientation);
                yMajor=centroid(2) - [-1 1]*(nuclearData{i,j}.MajorAxisLength/2)*sind(orientation);
                xMinor=centroid(1) - [-1 1]*(nuclearData{i,j}.MinorAxisLength/2)*sind(orientation);
                yMinor=centroid(2) - [-1 1]*(nuclearData{i,j}.MinorAxisLength/2)*cosd(orientation);
        
%                 colorRef=round(variable(i,j)*gradation)-round(minVariable*gradation);%define color
%     
%                 
%                 if colorRef==0 || isnan(colorRef)==1 %if value is absent, define to minimum color value
%         
%                     colorRef=1;
%         
%                 end
                if any(j==nuclei_subset)==1

                    for k=1:size(nuclearData{i,j}.PixelList,1)%Plot each pixel inside a given region with the appropriate color

                        plot(nuclearData{i,j}.PixelList(k,1),nuclearData{i,j}.PixelList(k,2),'.','Color','r','LineWidth',0.01);
        
                    end
                   
                end
                
                %plot convex hull defining outter boundary of segmented nucleus
                plot(nuclearData{i,j}.ConvexHull(:,1),nuclearData{i,j}.ConvexHull(:,2),'-k','LineWidth',2);
                
            end
    
        end
        
        %set plot formatting variables
        set(gca,'LineWidth',2.5,'FontSize',14);
        xlim([0 max(resolutions(:,1))])
        ylim([0 max(resolutions(:,2))])
        xlabel('x pixel')
        ylabel('y pixel')

        %colorbar;%show colorbar
        %caxis([minVariable maxVariable])
        title(strcat('t=',num2str(time(i)),'h'))
        
        %print figure into PDF file (change name as needed)
        %fname = sprintf('translational#1_concentration_heatmap_%d.pdf',time(i));
        print(h, '-dpdf', fname)
        %close(h)
        
end
