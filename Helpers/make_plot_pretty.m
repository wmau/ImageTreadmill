function [ hout] = make_plot_pretty(hin, varargin )
% hax_out = make_plot_pretty(hax_in, varargin )
%   Makes plots "pretty" (i.e. ready for minimal editing for a
%   manuscript/poster) by doing things like automatically making 
%   lineweights thicker for axes/lines/bar graphs/markers etc.
%
%   INPUTS
%   hin: handle to EITHER figure axes or graphics object(s). If axes it
%   will attempt to change all children of that axes (i.e. all lines, bars,
%   etc.)
%
%   PARAMETERS (optional, specific as ...'parameter1', value, 'parameter2, 
%               value2, etc.)
%
%   'type': graphics object type to include, e.g. 'line' or 'bar'
%
%   'linewidth': 2 is default
%
%   OUTPUTS
%   hout: modified hin

%% Parse Inputs
ip = inputParser;
ip.addRequired('hin',@(a) any(ishandle(a)));
ip.addParameter('type','all',@ischar)
ip.addParameter('linewidth', 2, @(a) isnumeric(a) && a > 0 && a <= 10);
ip.addParameter('fontsize', 14, @(a) isnumeric(a) && a > 0);
ip.addParameter('fonttype', 'SansSerif', @ischar);
ip.parse(hin,varargin{:})
type = ip.Results.type;
linewidth = ip.Results.linewidth;
fontsize = ip.Results.fontsize;
fonttype = ip.Results.fonttype;

%% Determine if axes only or just a single type of graphics object and break them apart
if length(hin) == 1 && isgraphics(hin,'axes')
    axes_flag = true;
    h_obj = get(hin,'Children');
    h_ax = hin;
    set(h_ax, 'Box', 'off') % Run top and right lines off on grid
    if ~strcmpi(type,'all') 
       h_obj = h_obj(isgraphics(h_obj,'type')); 
    end
else 
    axes_flag = false;
    h_obj = hin;
end

%% Modify axes and graphics objects
if axes_flag
    set(h_ax,'LineWidth',linewidth,'FontSize',fontsize,'FontName',fonttype);
end

if ~isempty(h_obj)
   for j = 1:length(h_obj)
       try
           set(h_obj(j),'LineWidth',linewidth);
       catch % ME - attempt to do fancy error catching - too much work
%            if strcmp(ME.message,  'There is no LineWidth property on the Group class.')
%                for j = 1:length(h_obj(j).Children
       end
   end
end


end
