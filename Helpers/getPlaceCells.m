function PCs = getPlaceCells(md,varargin)
%PCs = getPlaceCells(md, crit, varargin)
%
%   Gets all the cells with significant mutual spatial information. 
%
%   INPUTS
%       md: Session entry.
%
%       crit: p-value for consideration as a place cell. (default = 0.01)
%
%       optional...
%       nHits: min number of hits to be considered a place cell
%       (default = 10)
%
%       ratioHits: min ratio of passes through its field that a cell
%       must be active to be considered a place cell (default = 0.2);
%
%   OUTPUT
%       PCs: Indices of place cells. 
%
%

%% Parse inputs
ip = inputParser;
ip.addRequired('md', @isstruct);
ip.addOptional('crit', 0.01, @(a) a > 0 && a <= 1);
ip.addParameter('name_append','', @ischar);
ip.addParameter('nHits', 10, @(a) round(a,0) == a);
ip.addParameter('ratioHits', 0.2, @(a) a > 0 && a <= 1);
ip.KeepUnmatched = true; % pass along some of the variables above

ip.parse(md,varargin{:});
crit = ip.Results.crit;

name_append = ip.Results.name_append;
nHits = ip.Results.nHits;
ratioHits = ip.Results.ratioHits;

%% Get place cells.

    dirstr = ChangeDirectory(md.Animal, md.Date, md.Session);
    load(fullfile(dirstr,['Placefields' name_append '.mat']),'TMap_gauss','pval');
    load(fullfile(dirstr,['PlacefieldStats' name_append '.mat']),'PFnHits','bestPF','PFpcthits');
    load(fullfile(dirstr,['SpatialInfo' name_append '.mat']),'MI');

    
    [~,bestPF] = max(PFnHits,[],2);
    idx = sub2ind(size(PFnHits),1:size(PFnHits,1),bestPF');    
    PCs = find(pval<crit & MI'>0 & PFnHits(idx) > nHits & PFpcthits(idx) > ratioHits);

    
end