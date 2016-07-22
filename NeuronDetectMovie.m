function NeuronDetectMovie(md,movietype,clim,varargin)
%NeuronDetectMovie(MD,movietype,clim,varargin)
%
%   Writes a MP4 independent of TreadmillLog data. Takes an h5 file and
%   overlays detected neurons. 
%
%   INPUTS
%       MD: Session entry.
%
%       movietype: String, either 'd1','smoothed','dff', or 'slpdf'. 
%
%       clim: Color axis limits because brightness levels differ. 
%
%       varargin: String with accompanying input. Accepts...
%           'halfwindow' followed by a number (10, when comparing T2output
%           to d1 movies. Default = 0.
%           'alt_input' followed by variable name, indicating the variable
%           containing FT. Default = 'ProcOut.mat'.
%           'noi' (neurons of interest) followed by row indices of FT that
%           you want to highlight. Default = [1:size(FT,1)].
%

%% Setup.
    cd(md.Location); 
    
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x));
    p.addRequired('movietype',@(x) ischar(x));
    p.addRequired('clim',@(x) isnumeric(x)); 
    p.addParameter('halfwindow',10,@(x) isscalar(x)); 
    p.addParameter('alt_input','FinalOutput.mat',@(x) ischar(x)); 
    p.addParameter('noi',@(x) isnumeric(x));
    p.parse(md,movietype,clim,varargin{:});
    
    HalfWindow = p.Results.halfwindow;
    alt_input = p.Results.alt_input;
    
    load(fullfile(md.Location,alt_input),'FT');
    load(fullfile(md.Location,'ProcOut.mat'),'FT');
    try
        highlight = p.Results.noi;
    catch
        highlight = 1:size(FT,1);
    end
    
    %Type of processed movie.
    movietype = lower(movietype); 
    switch movietype
        case 'd1'
            h5file = fullfile(pwd,'D1Movie.h5');
            
            if ~exist(h5file,'file')
                cd(fullfile(md.Location,'MotCorrMovie-Objects'));
                motcorrh5 = dir('*.h5'); 
                
                TempSmoothMovie(fullfile(md.Location,'MotCorrMovie-Objects',motcorrh5.name),...
                    fullfile(md.Location,'SMovie.h5'),20);
                cd(md.Location);
                
                multiplier_use = DFDT_Movie('SMovie.h5','D1Movie.h5');
                if ~isempty(multiplier_use)
                    delete D1Movie.h5
                    multiplier_use = DFDT_Movie('SMovie.h5','D1Movie.h5',multiplier_use);
                    save multiplier.mat multiplier_use
                end
                delete SMovie.h5
            end

        case 'smoothed'
            cd('ICmovie_smoothed-Objects'); 
            h5file = fullfile(pwd,'Obj_1 - ICmovie_smoothed.h5'); 
            cd ..
        case 'dff'
            h5file = fullfile(pwd,'DFF.h5'); 
        case 'slpdf'
            h5file = fullfile(pwd,'SLPDF.h5'); 
    end
    
%% Write video. 
    %Set up for writing video. 
    imagingmovie = VideoWriter(['ImagingMovie_',movietype],'MPEG-4'); 
    imagingmovie.FrameRate = 20; 
    
    %For the for loop.
    info = h5info(h5file,'/Object'); 
    nFrames = info.Dataspace.Size(3); 
    
    %For plotting. 
    nNeurons = size(FT,1); 
    colors = rand(nNeurons,3); 
    outlines = cellfun(@bwboundaries,NeuronImage,'unif',0); 
    
    %Initialize writing. 
    open(imagingmovie); 
    ifigure = figure('Position',[260 240 560 420]); 
    
    %Write movie. 
    for i=1:nFrames
        %Get frame.
        frame = loadframe(h5file,i+HalfWindow,info); 

        %Find active neurons. 
        active = find(FT(:,i));
        
        %Display the imaging movie frame. 
        set(0,'currentfigure',ifigure); 
        imagesc(frame); caxis(clim); colormap gray; 
        
        %If there are active neurons according to FT, plot its outline. 
        if ~isempty(active) 
            hold on;
            for neuron=active'
                if ismember(neuron,highlight)
                    plot(outlines{neuron}{1}(:,2),outlines{neuron}{1}(:,1),...
                        'Color',colors(neuron,:),'linewidth',2);                   
                end
            end
            hold off;
        end
        
        %Get rid of axis marks. 
        axis equal; axis off; 

        %Write video to the brain imaging movie. 
        F = getframe(gcf); 
        writeVideo(imagingmovie,F); 
    end
                
    close(gcf); 
    close(imagingmovie); 
end