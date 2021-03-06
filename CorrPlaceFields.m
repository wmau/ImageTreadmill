function corrStats = CorrPlaceFields(ref,ssn,noi,varargin)
%corrStats = CorrPlaceFields(ref,ssn,noi,varargin)
%
%   Performs correlation on place field heatmaps. Also see CorrTrdmllTrace.

%% Parse inputs. 
    p = inputParser;
    p.addRequired('ref',@(x) isstruct(x));
    p.addRequired('ssn',@(x) isstruct(x)); 
    p.addRequired('noi',@(x) isnumeric(x)); 
    p.addParameter('pftype','placefieldsunsmoothed',@(x) ischar(x));
    p.addParameter('corrtype','pearson',@(x) ischar(x)); 
    p.addParameter('shuffle',false,@(x) islogical(x));
    
    p.parse(ref,ssn,noi,varargin{:});
    pftype = p.Results.pftype;
    corrtype = p.Results.corrtype;
    shuffle = p.Results.shuffle;

%% Match neurons. 
    ssns = [ref,ssn]; 
    DATA = CompileMultiSessionData(ssns,{pftype,'runoccmaps'});
   
    matchMat = msMatchCells(ssns,noi,true);       %Match cells.

    nNeurons = length(DATA.(pftype){1});
    corrStats = nan(nNeurons,2);
    
    maze1 = ~isnan(DATA.runoccmaps{1}(:));
    
%% Do correlations.
    for i = 1:size(matchMat,1)
        n1 = matchMat(i,1);
        
        if shuffle, n2 = matchMat(randsample(size(matchMat,1),1),2);
        else, n2 = matchMat(i,2); end
            
        %Placefields. 
        pf1 = DATA.(pftype){1}{n1}(:);
        pf2 = DATA.(pftype){2}{n2}(:);
        
        pf1 = pf1(maze1); %pf1(isnan(pf1)) = 0;
        pf2 = pf2(maze1); %pf2(isnan(pf2)) = 0;
        
        pf1 = pf1./max(pf1);
        pf2 = pf2./max(pf2);
        
        %Correlation.
        [corrStats(n1,1),corrStats(n1,2)] = corr(pf1,pf2,'type',corrtype,'rows','complete');
    
        if isnan(corrStats(n1,1)), corrStats(n1,1) = 0; corrStats(n1,2) = 1; end
    end
end