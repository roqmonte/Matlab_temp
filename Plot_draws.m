function hist = Plot_draws(results,info,print)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 17/Dec/2022
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Computes histograms for the draws of reduced form parameters
% and impulse response at time zero (impact matrices A0innv)
% Input:
%   results:
%   -.B_draws       : Draws from the reduced form parameters.
%   -.irf_full      : Impulse response functions from the model.
%   info:
%   -.names         : Vector with labels for variables.
%   -.shock_names   : Labels for shocks.
%   -.p             : Lag order.
%   -.determ_cte    : Contant term.
%   -.determ_trend  : Linear trend.
%   -.determ_trend2 : Quadratic term.
%   -.determ_exo    : Exogebous variable.
%   print           : Figure number
%
% Output:
%   No outputs
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Number of equaition and param per equation
neq  = size(results.B_draws,1);
npar = size(results.B_draws,2);

% Checking labels for shocks
nlab_sh = size(info.shock_names,2);
temp = neq-nlab_sh;
if temp > 0
    i = 1;
    for i0 = temp:-1:1
        info.shock_names(nlab_sh+i) = cellstr(strcat('Shcok','(',num2str(i),')'));
        i = i + 1;
    end
elseif temp < 0
    error('Wrong shock/label configuration') 
end

% Get labels for ylags
i = 1;
for i0 = 1:info.p
    for i1 = 1:neq
        temp = cellstr(strcat(char(info.names(i1)),'(',num2str(-i0),')'));
        labels_chart(i) = temp;
        i = i + 1;
    end
end

% Get labels for cte
if size(info.determ_cte,2) > 0
    labels_chart(i) = cellstr('Constant');
    i = i + 1;
end

% Get labels for linear trend and quadratic rend
if size(info.determ_trend,2) > 0
    labels_chart(i) = cellstr('Lin trend');
    i = i + 1;
end
if size(info.determ_trend2,2) > 0
    labels_chart(i) = cellstr('Squared trend');
    i = i + 1;
end

% Get labels for exo vars
if size(info.determ_exo,2) > 0
    for i1 = 1:size(info.determ_exo,2)
        temp = cellstr(strcat('exo var','(',num2str(i1),')'));
        labels_chart(i) = temp;
        i = i + 1;
    end
end
clear temp i i0 i1;

% Do histograms for draws of the reduced form parameters of the model
i = 1;
for i0 = 1:neq
    for i1 = 1:npar
        % Getting draws from patram bj
        bj_draws = squeeze(results.B_draws(i0,i1,:));
        
        % Define chart
        figure(print)
        subplot(neq,npar,i)
        
        % Histogram
        % Histogram
        if std(bj_draws) < 1e-4
            histogram(bj_draws.*0,50,'Normalization','probability')            
        else
            histogram(bj_draws,'Normalization','probability')
        end
        % Title for each column
        if i0 == 1
            title(labels_chart(i1),'FontSize',10);
        end        

        % Title for each row
        if i1 == 1
            ylab = char(info.names(i0));
            ylabel(ylab,'FontSize',12)
        end
        %set(gca,'FontSize',8);
        
        i = i + 1;
    end
end

% Do histograms for draws of the strectural matrix A0 (irf(0))
for i0 = 1:size(results.A0_draws,3)
    A0inv_draws(:,:,i0) = (results.A0_draws(:,:,i0))^(-1);
end
i = 1;
for i0 = 1:neq
    for i1 = 1:neq
        % Getting draws from A0(i0,i1)
        A0ij = squeeze(A0inv_draws(i0,i1,:));
        
        % Define chart
        figure(print+1)
        subplot(neq,neq,i)
        
        % Histogram
        if std(A0ij) < 1e-10
            histogram(A0ij,50,'Normalization','probability');
        else
            histogram(A0ij,'Normalization','probability');
        end
        % Title for each column
        if i0 == 1
            title(info.shock_names(i1),'FontSize',10);
        end        
        
        % Title for each row
        if i1 == 1
            ylab = char(info.names(i0));
            ylabel(ylab,'FontSize',12)
        end
        
        
        i = i + 1;
    end
end
hist = [];
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%