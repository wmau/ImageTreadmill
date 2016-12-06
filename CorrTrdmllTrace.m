function corrStats = CorrTrdmllTrace(ref,ssn,noi,varargin)
%corrStats = CorrTrdmllTrace(ref,ssn,noi,corrtype)
%
%   Performs correlation on time field tuning curves. Also see
%   CorrPlaceFields. 

%% Parse inputs. 
    p = inputParser;
    p.addRequired('ref',@(x) isstruct(x));
    p.addRequired('ssn',@(x) isstruct(x)); 
    p.addRequired('noi',@(x) isnumeric(x)); 
    p.addParameter('corrtype','pearson',@(x) ischar(x)); 
    p.addParameter('tracetype','curves',@(x) ischar(x)); 
    
    p.parse(ref,ssn,noi,varargin{:});
    corrtype = p.Results.corrtype;
    tracetype = p.Results.tracetype; 
    
%% Get mapped neurons.
    ssns = [ref,ssn];   
    DATA = CompileMultiSessionData(ssns,{tracetype});
   
    mapMD = getMapMD(ref);
    matchMat = msMatchCells(mapMD,ssns,noi,true);
    
    nNeurons = length(DATA.curves{1}.tuning);
    corrStats = nan(nNeurons,2);
    
%% Do correlations.
    for i = 1:size(matchMat,1)
        n1 = matchMat(i,1);
        n2 = matchMat(i,2);
        
        switch tracetype
            case 'curves'
                tf1 = DATA.curves{1}.tuning{n1}';
                tf2 = DATA.curves{2}.tuning{n2}';
            case {'rawtrdmll','difftrdmll','tracetrdmll'}
                tf1 = mean(DATA.(tracetype){1}(:,:,n1))';
                tf2 = mean(DATA.(tracetype){2}(:,:,n2))';
        end
        
        [corrStats(n1,1),corrStats(n1,2)] = corr(tf1,tf2,'type',corrtype);
    end
end