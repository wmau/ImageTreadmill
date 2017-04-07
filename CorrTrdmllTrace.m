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
    p.addParameter('shuffle',false,@(x) islogical(x));
    
    p.parse(ref,ssn,noi,varargin{:});
    corrtype = p.Results.corrtype;
    tracetype = p.Results.tracetype; 
    shuffle = p.Results.shuffle;
    
%% Get mapped neurons.
    ssns = [ref,ssn];   
    DATA = CompileMultiSessionData(ssns,{tracetype});
   
    mapMD = getMapMD(ref);
    matchMat = msMatchCells(mapMD,ssns,noi,true);
    
    try nNeurons = length(DATA.curves{1}.tuning);
    catch, nNeurons = size(DATA.(tracetype){1},3); end
    corrStats = nan(nNeurons,2);
    
%% Do correlations.
    for i = 1:size(matchMat,1)
        n1 = matchMat(i,1);
        
        if shuffle, n2 = matchMat(randsample(size(matchMat,1),1),2);
        else, n2 = matchMat(i,2); end
        
        switch tracetype
            case 'curves'
                tf1 = DATA.curves{1}.tuning{n1}';
                tf2 = DATA.curves{2}.tuning{n2}';
            case {'rawtrdmll','dfdttrdmll','lptrdmll'}
                tf1 = mean(DATA.(tracetype){1}(:,:,n1))';
                tf2 = mean(DATA.(tracetype){2}(:,:,n2))';
        end
        
%         tf1 = tf1./max(tf1);
%         tf2 = tf2./max(tf2);
        [corrStats(n1,1),corrStats(n1,2)] = corr(tf1,tf2,'type',corrtype);
        
        if isnan(corrStats(n1,1)), corrStats(n1,1) = 0; corrStats(n1,2) = 1; end
    end
end