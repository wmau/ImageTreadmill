function [all,notime,nocells,automated] = SourceSinkGLM(tbl,tracetype)
%[all,notime,nocells,automated] = SourceSinkGLM(tbl)
%
%   Fits GLMs to the table in multiple ways. 
%       -Using all the regressors. 
%       -Excluding time.
%       -Excluding cells. 
%       -Stepwise regression with AIC as the parameter addition/removal
%       criterion.
%
%   INPUT
%       tbl: Table of predictor variables with variable names in the form
%       'n#' where # is the neuron # in FT. The first column is time and
%       the last column is the response variable. 
%
%   OUTPUTS
%       all: GLM containing all the predictors including time. 
%
%       notime: GLM containing all neural predictors, excluding time. 
%
%       nocells: GLM containing only time as a predictor. 
%
%       automated: 

%% GLM fits
    if strcmp(tracetype,'FT'), dtype = 'binomial'; linkfx = 'logit'; 
    elseif strcmp(tracetype,'rawtrace'), dtype = 'normal'; linkfx = 'identity'; end
    all = fitglm(tbl,'distribution',dtype,'link',linkfx);
    notime = fitglm(tbl(:,2:end),'distribution',dtype,'link',linkfx);
    nocells = fitglm(tbl(:,[1,end]),'distribution',dtype,'link',linkfx);
    automated = stepwiseglm(tbl,'linear',...
         'upper','linear',...
         'distribution',dtype,...
         'link',linkfx,...
         'Criterion','AIC');
end