function d = cohensD(x,y)
%d = cohensD(x,y)
%
%   Name says it all. 
%
%   INPUTS
%       X & Y: two vectors that you want to Cohen's d. 
%
%   OUTPUT
%       d: Obvious.
%

%% Do the Cohen's d.
    M1 = mean(x);
    M2 = mean(y);
    SDpooled = sqrt((std(x)^2 + std(y)^2)/2);
    
    d = (M1-M2)/SDpooled;
    
end