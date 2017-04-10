function Alt = postrials_treadmill(md,varargin)
%function data = postrials(x,y,plot_each_trial,...)
%   
%   This function takes mouse position data and sorts them into trial
%   numbers, left/right, and correct/incorrect.
%
%   INPUTS: 
%       X & Y: Position vectors after passing through
%       PreProcessMousePosition.
%       
%       plot_each_trial: Logical to determine whether or not you want the
%       function to plot the XY position of the mouse for each trial. 
%       
%   OUTPUTS:
%       DATA: a struct with these fields:
%           frame = Frame number.
%           trial = Trial number. 
%           choice = Left (1) vs. right (2).
%           alt = Correct (1) vs. incorrect (0). 
%           x = X position.
%           y = Y position.
%           section = section number. Refer to getsection.m. 
%           goal = goal location. 1 = left. 2 = right. 0 = not goal
%           location. 
%           summary = Summary of trials. The first column is trial number
%           followed by left/right and correct/incorrect in the same format
%           as above. 
%
%   TIP: To find frames for a particular trial of interest, you can do:
%       data.frames(data.trial == TRIAL_OF_INTEREST).
%
%%
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x));
    p.addParameter('plotEachTrial',false,@(x) islogical(x));
    p.addParameter('suppressOutput',true,@(x) islogical(x)); 
    
    p.parse(md,varargin{:});
    
    plotEachTrial = p.Results.plotEachTrial; 
    suppressOutput = p.Results.suppressOutput;
    
%% Get XY coordinates.
    cd(md.Location);
    load('Pos_align.mat','x_adj_cm','y_adj_cm');
    x = x_adj_cm; y = y_adj_cm; 
    
    load('TimeCells.mat','TodayTreadmillLog');
    direction = TodayTreadmillLog.direction; 
    
    [~,folder] = fileparts(md.Location); 
    blocked = ~isempty(strfind(folder,'blocked'));  

%% Label position data with section numbers. 
    bounds = sections_treadmill(x,y,direction,0);
    [sect,rot_x,rot_y] = getsection_treadmill(x,y,bounds);
    
%% Define important section numbers. 
    %Define sequences of section numbers that correspond to left or right
    %trials. **CURRENTLY OBSOLETE**
    %left = 1:6; 
    %right = [1,2,3,7,8,9]; 
    
    %Define opposite arms. Important for catching miscategorization of 
    %trials by this script. 
    left = 5; right = 8; 
    
    %Define return arms. **CURRENTLY OBSOLETE**
    %return_arm = [6,9]; 
    
    %Define when the mouse is at the start position. 
    base = find(sect==1);
    
%% Label trials. 
    %Define first trial. When does the mouse first enter the starting
    %location? 
    start = min(find(sect==1)); 

    %Preallocate.
    epochs = start; 
    mouse_running = 1;
    
    %For each lap. 
    for this_trial = 1:200
        
        try     %Try sorting a trial. 
        
        %Index for next trial. 
        next = this_trial+1;    
        epochs(next) = epochs(this_trial)+1;    %First glance: trials are at least 1 frame long. 

        %As long as the criteria are not satisfied (see below)...
        while mouse_running
            
            epochs(next) = epochs(next) + 1;    %...keep adding frames.

            %To break the while loop, must satisfy the criteria that the
            %mouse entered the left/right arm from the left/right approach
            %area. 
            if (sect(epochs(next)) == left && sect(epochs(next)-1) == left-1) || ...    
                    (sect(epochs(next)) == right && sect(epochs(next)-1) == right-1)
                
                %Label this trial as left/right.
                if sect(epochs(next)) == left
                    trialtype(this_trial) = 1; 
                elseif sect(epochs(next)) == right
                    trialtype(this_trial) = 2; 
                end
                
                %Once mouse enters left/right arm, reach for next instance
                %where he is in the base. The duration since the while loop
                %started up until this timepoint is now one lap. 
                try
                    epochs(next) = min(base(base > epochs(next))); 
                catch
                    epochs(next) = length(rot_x); 
                end
             
                %Plot laps. 
                if plotEachTrial
                    figure(this_trial);
                    plot(x(epochs(this_trial):epochs(next)), y(epochs(this_trial):epochs(next))); 
                    xlim([min(x) max(x)]); ylim([min(y) max(y)]); 
                    title(['Trial ', num2str(this_trial)], 'fontsize', 12); 
                end

                %Notify user of possible errors in the trial sorting script. This
                %catches when the mouse appears on both maze arms in what the
                %script believed to be a single trial. 
                if ((ismember(left, sect(epochs(this_trial):epochs(next))) && trialtype(this_trial) == 2) || ...
                        (ismember(right, sect(epochs(this_trial):epochs(next))) && trialtype(this_trial) == 1)) && ...
                        ~suppressOutput
                    disp(['Warning: This epoch may contain more than one trial: Trial ', num2str(this_trial)]);
                end       
                
                break;
            end
        end
        
        %When postrials can no longer successfully sort a trial, stop the
        %loop. 
        catch
            numtrials = this_trial - 1;             
            if length(trialtype) > numtrials    %Sometimes, the mouse leaves the maze on a left/right arm so trialtype gets an extra entry.
                trialtype(this_trial) = [];     %In that case, remove it. 
            end
            break; 
        end
        
    end
    
%% Build up the struct. 
    Alt.frames = 1:length(x);          %Frames.
    
    %Vector containing correct vs. error using a trick: take the difference
    %between consecutive trial types (left (1) vs. right (2)) such that
    %errors corresponding to consecutive double visits will be zero (one
    %minus one or two minus two). Take the absolute value and all other
    %correct visits will be 1. The first choice (reward on both arms) is
    %always correct. 
    if blocked 
        alt = ones(1,length(trialtype));
    else
        alt = [1 abs(diff(trialtype))]; 
    end
    
    %Trial numbers, trial type (left vs. right), and correct vs. incorrect. 
    for this_trial = 1:numtrials
        next = this_trial+1; 

        Alt.trial(epochs(this_trial):epochs(next)) = this_trial; 
        Alt.choice(epochs(this_trial):epochs(next)) = trialtype(this_trial); 
        Alt.alt(epochs(this_trial):epochs(next)) = alt(this_trial); 
    end
    
    %This script may not cover the entire session. If that's the case, pad
    %the rest of it with 0s or NaNs.
    if length(Alt.trial) < length(x)
        Alt.trial(end:length(x)) = 0; 
        Alt.choice(end:length(x)) = 0;
        Alt.alt(end:length(x)) = NaN;      %NaNs so that the last uncaptured trials don't register as errors. 
    end
    
    %Mouse position. 
    Alt.x = rot_x;                  %X position.
    Alt.y = rot_y;                  %Y position. 
    Alt.section = sect';            %Section number. Refer to getsection.m.
    
    %Summary. 
    Alt.summary = [(1:numtrials)', trialtype', alt']; 
    
    %Save. 
    save Alternation Alt;
end