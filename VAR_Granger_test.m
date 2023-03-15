function results = VAR_Granger_test(y1,x1,trend,exo,max_lag,Opt)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 14/Mar/2022
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Granger test using VAR functions. The optimizes the lag
% selection of the test and performs likelihood ratio test. Constant term
% added by default.
% Inputs:
%   y1          : Variable y
%   x2          : Variable x
%   trend       : (0) No trend; (1) linear trend; (2) linear and quadratic trend.
%   exo         : Exogenous variables beside cte.
%   max_lag     : Max lag for optimizer
%   Opt         : Selection criteria:(1) AIC,(2) HQC default , and (3) BIC.
%
% Outputs:
%   results:
%   -.GT_x_y    : Test H0: xi -> yi. Gt test and Pvalue
%   -.GT_y_x    : Test H0: yi -> xi. Gt test and Pvalue
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Checking inputs
% If no exo var, define cte term
if exist('exo','var') == 0 || isempty(exo)
    exo = [];
end
% Checking trend
if exist('trend','var') == 0 || isempty(trend)
    trend = 0;
end
% If no Opt, use HQC
if exist('Opt','var') == 0 || isempty(Opt)
    Opt = 1;
end

% Basic set up for model
cte    = 1;
info.p = max_lag;
info.px= 0;
% Setting dates.
info.dates_ini  = [2010,1,1];       % First observation of the data.
% Settings for impulse responses.
% Options for plots.
info.names      = {'Y1 var' 'X2 var'};

% Preparing data for estimation.
data_ini = [y1 x1];
[data,determ] = data_make(data_ini,info,cte,trend,exo);

% Testing lag lenght of the model.
info.alpha   = 0.05;
info.max_lag = max_lag;
lag_res      = VAR_TestLagLength(data,info,determ);
% Selection criteria
if Opt== 1
    [IC,id] =min(lag_res.aic');
elseif Opt== 2
    [IC,id] =min(lag_res.hqc');
elseif Opt==3
    [IC,id] =min(lag_res.sic');
end
% Fixing if no lag is selected
if id == 1
    p = id;
else
    p = id-1;
end
clear id IC;

% Getting data ready
yfull = LagN(data.all(:,1),p);
xfull = LagN(data.all(:,2),p);
% Test H0: xi -> yi
[GT_x_y,Pv_x_y] = granger(xfull,yfull,determ(1+p:end,:));
% Test H0: yi -> xi
[GT_y_x,Pv_y_x] = granger(yfull,xfull,determ(1+p:end,:));

% Saving results
results.GT_x_y = [GT_x_y,Pv_x_y];
results.GT_y_x = [GT_y_x,Pv_y_x];

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Functions
% Preparing data for VAR model estimation and lag order selection
function [data,determ] = data_make(data_ini,info,cte,trend,exo)
% Description: Generates data for VAR model
% Input:
%   data_ini        : Matrix with all data.
%   info:
%   -.px            : Lag order exo variables.
%   cte             : (0) No constant; (1) Constant (default).
%   trend           : (0) No trend (default); (1) linear trend; (2) quadratic trend.
%   exo             : Exogenous variables of the model.
%
% Output:
%   data            : Data for model
%   determ          : Constant term, trends and exo variables for VAR.
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Checking exo variables
if exist('exo','var') == 0
    exo = [];
    px = 0;
else
    px = info.px*(1-isempty(exo));
end

% Results.
data.endo  = data_ini(1+px:end,:);
data.all   = data_ini(1+px:end,:);

% Building exo variables.
% Constant term
if cte == 0
    info.determ_cte = [];
elseif cte == 1
    info.determ_cte = ones(size(data_ini,1),1);    
else
   error('Check cte variable'); 
end
% Trends
if trend == 0
    info.determ_trend = [];
    info.determ_trend2= [];
elseif trend == 1
    info.determ_trend = (1:size(data_ini,1))';
    info.determ_trend2= [];
elseif trend == 2
    info.determ_trend = (1:size(data_ini,1))';
    info.determ_trend2= ((1:size(data_ini,1)).^(2))';
else
   error('Check trend variable'); 
end
% Exo variables.
if size(exo,2) == 0 
    info.determ_exo = [];
elseif size(exo,2) > 0
    exo2 = LagN(exo,px);
    info.determ_exo = exo2;
end
% Deterministic variables.
determ = [info.determ_cte(1+px:end) info.determ_trend(1+px:end) info.determ_trend2(1+px:end) info.determ_exo];

% Do GT test for H0: yi -> xi
function [Ftest,Pval] = granger(yfull,xfull,determ)
% Restricted model for Xi
x    = xfull(:,1);
xlags= xfull(:,2:end);
kx   = size(xlags,2);
% OLS Estimation
mod_res = OLSest(x,[xlags determ]);
SSR_res = mod_res.SSR;

% Unretsrcited model
y    = yfull(:,1);
ylags= yfull(:,2:end);
ky   = size(ylags,2);
% OLS Estimation
mod_unres = OLSest(x,[xlags ylags determ]);
SSR_unres = mod_unres.SSR;

% Computing stat
q = mod_unres.k - mod_res.k;
T =size(x,1);
K = mod_unres.k;

Ftest = ((SSR_res - SSR_unres) / SSR_unres)* ((T-K)/q);
Pval  = 1 - fcdf(Ftest,q,(T-K));

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%