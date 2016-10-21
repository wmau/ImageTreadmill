function [t,i] = getTimePeak(MD)
%[t,i] = getTimePeak(MD)
%
%   Loads tuning curves and searches for all the peaks. 
%
%   INPUT
%       MD entry.
%
%   OUTPUTS
%       t: Time peak. 
%
%       i: Index peak. 
%

%%
    load(fullfile(MD.Location,'TimeCells.mat'),'curves','T');
    
    [~,i] = cellfun(@max,curves.tuning);
    t = i/4; 
end