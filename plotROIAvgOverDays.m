function plotROIAvgOverDays(referenceMD,otherMDs,neurons,varargin)
%plotROIAvgOverDays(referenceMD,otherMDs,neurons,varargin)
%
%   Plots average fluorescence during transients of specified cells across
%   days. List the neuron numbers as they appear in referenceMD. The
%   function will draw a 20-pixel window around the cell centered on its
%   centroid and then use the same window on subsequent sessions to plot
%   cell ROIs. 
%
%   INPUTS
%       referenceMD: session whose neuron indices you want to use. 
%
%       otherMDs: other sessions to visualize specified neurons. 
%
%       neurons: neurons you want to plot. 
%
%       varargin: 
%           plotOrientation, 'vertical' or 'horizontal': orientation of
%           subplots. 
%       
%           scalebar: logical, plots a 10 micron scale bar. 
%
%   OUTPUT
%       Plot of ROIs over days. 
%

%% Parse inputs. 
    p = inputParser; 
    p.addRequired('referenceMD',@(x) isstruct(x)); 
    p.addRequired('MDs',@(x) isstruct(x));
    p.addRequired('neurons',@(x) isnumeric(x)); 
    p.addParameter('plotOrientation','horizontal',@(x) ischar(x));
    p.addParameter('scalebar',true,@(x) islogical(x)); 
    
    p.parse(referenceMD,otherMDs,neurons,varargin{:});
    plotOrientation = p.Results.plotOrientation;
    scalebar = p.Results.scalebar; 
    
    nNeurons = length(neurons);
    
    cd(referenceMD.Location);
    load('FinalOutput.mat','NeuronImage'); 
    dims = size(NeuronImage{1});                %Dimensions of FOV. 
    
%% Collect initial information.
    %Combine session entries. Also get session that contains cell mapping. 
    MDs = [referenceMD otherMDs];
    nSessions = length(MDs); 
   
    %Sort the dates so that we plot in chronological order. 
    dates = {MDs.Date}; 
    [sortedDates,chronOrder] = sort(dates);
    
    %Replace the underscores with dashes. 
    for s=1:nSessions
        sortedDates{s}(3:3:6) = '-';
    end
    
%% Collect neuron indices on different sessions and their registration info.
    matchedCells = msMatchCells(MDs,neurons,false);
    
    %Get registration information.
    for s=2:nSessions
        regInfo(s) = image_registerX(referenceMD.Animal,referenceMD.Date,...
            referenceMD.Session,otherMDs(s-1).Date,otherMDs(s-1).Session);
    end
    
%% Build the ROIs. 
    %Preallocate. 
    [centroids,allNeuronAvg] = deal(cell(nSessions,1)); 
    [xWindow,yWindow] = deal(cell(nNeurons,1)); 
    for s=1:nSessions
        cd(MDs(s).Location);
        load('FinalOutput.mat','NeuronImage','NeuronAvg'); 
        
        %Make the average mask. 
        allNeuronAvg{s} = MakeAvgROI(NeuronImage,NeuronAvg);  
        
        %Get centroids.
        props = cellfun(@(x) regionprops(x,'Centroid'),NeuronImage);
        temp = extractfield(props,'Centroid'); 
        centroids{s} = [temp(1:2:end)', temp(2:2:end)'];
        
        %For sessions that aren't the base session, transform the masks so
        %they align. 
        if s~=1
            allNeuronAvg{s} = cellfun(@(x) imwarp_quick(x,regInfo(s)),...
                allNeuronAvg{s},'UniformOutput',false); 
        end
        
        %Make a window around the neuron.        
        if s==1 
            for n=1:nNeurons
                %Neuron number. 
                thisNeuron = matchedCells(n,s);  
                
                %Columns in the image. 
                xLowerBound = round(centroids{s}(thisNeuron,1)-10); 
                xUpperBound = xLowerBound + 20;     
                
                %To handle cases on the edge. 
                if xLowerBound < 1, xLowerBound = 1; end                %Ran out of space on the left?
                if xUpperBound > dims(2), xUpperBound = dims(2); end    %Ran out of space on the right?
            
                %Make the window. 
                xWindow{n} = xLowerBound:xUpperBound;

                %Rows in the image. 
                yLowerBound = round(centroids{s}(thisNeuron,2)-10);     %Bottommost boundary.
                yUpperBound = yLowerBound + 20;                         %Topmost boundary. 

                %To handle cases on the edge. 
                if yLowerBound < 1, yLowerBound = 1; end                %Ran out of space at the bottom? 
                if yUpperBound > dims(1), yUpperBound = dims(1); end    %Ran out of space at the top?

                %Make the window. 
                yWindow{n} = yLowerBound:yUpperBound; 
            end   
        end      
    end    
    
    %For each neuron specified, draw a window around it centered on the
    %centroid of the reference session. 
    for n=1:nNeurons       
        for s=1:nSessions
            %Neuron index. 
            thisNeuron = matchedCells(n,s); 
            
            %Window. 
            if thisNeuron > 0 
                allNeuronAvg{s}{thisNeuron} = ...
                    allNeuronAvg{s}{thisNeuron}(yWindow{n},xWindow{n});
            end
        end
    end
    
%% Plot. 
    keepgoing = true;
    n = 1;
    while keepgoing        
        f = figure(1216); hold on;
        plotNumber = 1;
        for s=chronOrder
            %Neuron index. 
            thisNeuron = matchedCells(n,s);
            
            %Get ROI. 
            if thisNeuron > 0, ROI = allNeuronAvg{s}{thisNeuron}; end

            %Plot either vertically or horizontally. 
            switch plotOrientation
                case 'horizontal', subplot(1,nSessions,plotNumber);
                case 'vertical', subplot(nSessions,1,plotNumber); 
            end
            
            %Plot the ROI. 
            imagesc(flipud(ROI));
            if scalebar
                line([0 10*1.10],[2 2],'linewidth',5,'color','w');
            end
            colormap gray;
            axis equal; axis off; 
            set(gca,'ydir','reverse'); 
            title({['Cell #',num2str(thisNeuron)],...
                sortedDates{s}});
            plotNumber = plotNumber + 1;
        end
        
        %Scroll through neurons. 
        [keepgoing,n] = scroll(n,nNeurons,f);
        close all;
    end
    
end