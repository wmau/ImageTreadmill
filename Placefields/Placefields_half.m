function [ ] = Placefields_half( MD, calc_mode, exclude_frames, name_append, varargin )
%Placefields_half(MD, mode, exclude_frames, name_append, varargin)
%
%   Wrapper function that runs Placefields on 1st v 2nd half or odd v even
%   minute data.  Must run Placefields on full data set first!
%
%   INPUTS
%       MD: session to analyze. 
%
%       mode: 'half' or 'oddeven'.
%
%       exclude_frames: bad frames to exclude. Should correspond
%       to bad frames in imaging movie / PSAbool from FinalOutput.mat
%
%       name_append: name to append to save files (optional).
%
%       optional...
%       See Placefields.

%% Parse inputs.

    [dirstr, MD] = ChangeDirectory(MD.Animal, MD.Date, MD.Session, 0); % Change Directory and fill in partial MD if used
    
    ip = inputParser;
    ip.addRequired('MD',@(x) isstruct(x)); 
    ip.addRequired('calc_mode',@(a) ischar(a) && (strcmpi(a,'half') || strcmpi(a,'oddeven')));
    ip.addRequired('exclude_frames', @(x) isnumeric(x)); 
    ip.addRequired('name_append', @ischar);
    ip.KeepUnmatched = true;
    % Note that other varargins will be checked in Placefields function
    
    ip.parse(MD, calc_mode, exclude_frames, name_append, varargin{:});

    %% Calculate indices for each half
    load(fullfile(dirstr,'FinalOutput.mat'),'PSAbool');
    load(fullfile(dirstr,'Pos.mat'),'time_interp');
    Flength = min([size(PSAbool,2), length(time_interp)]);
    
    ind_use_half{1} = false(1,Flength);
    ind_use_half{2} = false(1,Flength);
    switch calc_mode
        
        case 'half' % Get indices to include for each half
            half = round(Flength/2);
            ind_use_half{1}(1:half) = true; 
            ind_use_half{2}(half+1:Flength) = true; 
        case 'oddeven'
            SR = 20; %fps for imaging acquisition
            epoch_chunk = SR*60; % # frames/minute
            
            last_frame = 0;
            epoch_add = 'odd';
            even_frames = [];
            odd_frames = [];
            while last_frame < Flength
               if strcmpi(epoch_add,'odd')
                   [odd_frames, last_frame] = build_epoch(odd_frames, epoch_chunk, last_frame, Flength);
                   epoch_add = 'even';
               elseif strcmpi(epoch_add,'even')
                   [even_frames, last_frame] = build_epoch(even_frames, epoch_chunk, last_frame, Flength);
                   epoch_add = 'odd';
               end
            end
            
            ind_use_half{1}(odd_frames) = true;
            ind_use_half{2}(even_frames) = true;
        otherwise
    end
    
    % Incorporate excluded frames
    include_frames_logical = true(1,Flength);
    include_frames_logical(exclude_frames) = false;
    
    exclude_half = cell(1,2);
    for j = 1:2
        ind_use_half{j} = ind_use_half{j} & include_frames_logical;
        exclude_half{j} = find(~ind_use_half{j});
    end

    %% Run Placefields for each half
    Placefields_halves = cell(1,2);
    for j = 1:2
        disp(['Calculating PFs for half # ' num2str(j)])
        Placefields(MD, 'name_append', [name_append '_' calc_mode num2str(j)], ...
            'exclude_frames', exclude_half{j}, varargin{:});
        PFfilename = fullfile(dirstr,['Placefields' name_append '_' calc_mode num2str(j) '.mat']);
        Placefields_halves{j} = load(PFfilename);
        delete(PFfilename);
    end
    
    save(fullfile(dirstr,['Placefields' name_append '_' calc_mode '.mat']),'Placefields_halves');
end

%% Build odd and even epoch frames
function [epoch_frames_out, last_frame_out] = build_epoch(epoch_frames_in, chunk, last_frame, num_frames)

epoch_frames_out = [epoch_frames_in, (last_frame+1):1:min([(last_frame+chunk), num_frames])];
last_frame_out = epoch_frames_out(end);

end
