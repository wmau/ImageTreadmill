function [keepgoing,i] = scroll(i,N,f)
%[keepgoing,i] = scroll(i,N,f)
%   
%   Function for scrolling through plots of arbitrary things. Press left
%   key to go back and right key to go forward. Esc to exit. 
%
%   INPUTS
%       i: Index of whatever you're plottnig. 
%
%       N: Total number of things you're plotting. 
%
%       f: Figure handle.
%
%   OUTPUTS
%       keepgoing: For a while loop. 0 to stop, 1 to keep going. 
%
%       i: New index after adding or subtracting one. 
%

%% Main function. 
    f; 
    
    [~,~,key] = ginput(1); 

    if key == 29 && i < N
        i = i + 1; 
        keepgoing = 1; 
    elseif key == 28 && i ~= 1
        i = i - 1; 
        keepgoing = 1; 
    elseif key == 27
        keepgoing = 0; 
        close(f); 
    else
        keepgoing = 1;
    end
    
end