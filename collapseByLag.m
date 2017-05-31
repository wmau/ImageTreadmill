function [m,sem,diags] = collapseByLag(R)
%[m,sem] = collapseByLag(R)
%
%   Computes the mean of the diagonals in the correlation coefficient
%   matrix R. 
%
%   INPUT
%       R: Correlation coefficient matrix. This is usually 5x5xZ where Z is
%       the number of instances (could be trials, sessions, or animals). 
%

%% Set up.
    [nRows,nCols,nInstances] = size(R);     
    nLags = min([nRows nCols]);
    diags = cell(nLags,1);              %Contains the diagonal elements of R slices. 

%% Get the means.
    %An "instance" here can refer to an animal, session, or trial. They are
    %iterations across which you will take the mean. This is always going
    %to be the third dimension in your R matrix. 
    for i=1:nInstances
        %Get a slice of the R matrix. 
        Rinstance = R(:,:,i);
        
        %For each time lag (usually 5)...
        for l=1:nLags
            %Get the row and column indices then covert them into linear. 
            row = 1:nLags-l+1;
            col = l:nLags;          %Look carefully -- that's an L.
            inds = sub2ind([nLags,nLags],row,col);
            
            %Get the diagonal of the R slice. 
            if isnumeric(Rinstance)
                diagOfR = Rinstance(inds);
            elseif iscell(Rinstance)
                diagOfR = cell2mat(Rinstance(inds)')';
            end
            
            %Append it to the cell array. 
            diags{l} = [diags{l} diagOfR];
        end
    end
    
    %Take the mean and SEM of the diagonals. 
    m = cellfun(@nanmean, diags);
    sem = cellfun(@nanstd, diags)./sqrt(cellfun('length', diags));
end
    