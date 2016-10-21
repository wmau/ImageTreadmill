function corrStats = CorrPlaceFields(ref,ssn,noi,varargin)
%corrStats = CorrPlaceFields(ref,ssn,noi,varargin)
%
%   

%% Parse inputs. 
    p = inputParser;
    p.addRequired('ref',@(x) isstruct(x));
    p.addRequired('ssn',@(x) isstruct(x)); 
    p.addRequired('noi',@(x) isnumeric(x)); 
    p.addParameter('pftype','placefieldsunsmoothed',@(x) ischar(x));
    p.addParameter('corrtype','pearson',@(x) ischar(x)); 
    
    p.parse(ref,ssn,noi,varargin{:});
    pftype = p.Results.pftype;
    corrtype = p.Results.corrtype;

%% Match neurons. 
    ssns = [ref,ssn]; 
    DATA = CompileMultiSessionData(ssns,{pftype});
   
    mapMD = getMapMD(ref);                      %Find reference map.
    matchMat = msMatchCells(mapMD,ssns,noi);    %Match cells.
    matchMat(matchMat(:,2)==0,:) = [];          %Get rid of unmatched neurons.
    matchMat(isnan(matchMat(:,2)),:) = [];      %Get rid of unmatched neurons.
    
    nNeurons = length(DATA.(pftype){1});
    corrStats = nan(nNeurons,2);
    
%% Do correlations.
    for n1 = 1:size(matchMat,1)
        n2 = matchMat(n1,2);
        
        %Placefields. 
        pf1 = DATA.(pftype){n1};
        pf2 = DATA.(pftype){n2};
        
        %Correlation.
        [corrStats(n1,1),corrStats(n2,2)] = corr(pf1,pf2,'type',corrtype);
    end
end