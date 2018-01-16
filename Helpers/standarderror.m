function SEM = standarderror(data)
%SEM = standarderror(data)
%
%   Calculates the standarde error of the mean. 
%

%% Compute standard error.
    SEM = nanstd(data)/sqrt(length(data(~isnan(data))));
end
    