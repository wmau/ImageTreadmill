function SEM = standarderror(data)
%SEM = standarderror(data)
%
%   Calculates the standarde error of the mean. 
%

%% Compute standard error.
    SEM = std(data)/sqrt(length(data));
end
    