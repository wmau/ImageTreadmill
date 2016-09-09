function output = sortbyDate(dates,input)
%output = sortDate(dates,input)
%   
%   Sorts an input chrononlogically based on corresponding dates. Can sort
%   vectors and matrices (by column).
%
%   INPUTS
%       dates: Cell array of dates in the format 'mm_dd_yyyy'.
%
%       input: Should align with dates. If it is a vector, must be same
%       size as dates. If it is a matrix, each column must corresponding to
%       a date. 
%

%% Sorting.  
    [~,order] = sort(datenum(dates,'mm_dd_yyyy'));
    if ndims(input)==1
        assert(length(dates)==length(input),...
            'Number of dates and size of input vector do not match!');
        output = input(order);
    elseif ndims(input)==2
        assert(length(dates)==size(input,2),...
            'Number of dates and number of columns of input matrix do not match!');
        output = input(:,order); 
    end
    
end