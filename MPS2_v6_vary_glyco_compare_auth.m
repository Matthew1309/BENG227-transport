% MPS2_improvement.m
% Source: https://www.sciencedirect.com/science/article/pii/S0022519321003027
% Purpose:
% Take the author's simple model of transport across the endothelial
% barrier and improve upon it with the details of a glycocalyx layer.
% Authors: Matthew Kozubov, Cal McKinney
% This version almost has it mimicking the author's data, but not entirely

%% Read in the author data
author_fig16 = readmatrix("phenomenological-based-model-adipose-tissue-glucose-levels.csv");
author_fig16 = author_fig16(2:end,:);
author_fig16(:,1) = author_fig16(:,1) - author_fig16(1,1);
plot(author_fig16(:,1), author_fig16(:,2));

%% Making a function that can use our optimized params
close all
%       1.4064
theta = [1.4064, 1.68*10^3];

y0 = [5.73e-3,8e-4,4e-3,4.7e-3]; %femto-mol/um^3
tspan = linspace(0,2450);
[t,Code45] = ode45( @(t,y)matcal_system(t,y,theta), tspan, y0);
figure(1)
hold on
h1 = plot(t,Code45(:,1));
h3 = plot(t,Code45(:,2));
h4 = plot(t,Code45(:,3));
h5 = plot(t,Code45(:,4));
% h2 = plot(t, yvals2, 'o', 'MarkerEdgeColor', 'magenta', 'LineStyle', 'none');
% legend([h1,h3,h4,h5,h2], ["Artery", "Glyco", "EC", "Adipose", "Auth data"]);
hold off

figure(2)
hold on
yvals2 = 10^-3 * interp1(author_fig16(:,1),author_fig16(:,2),tspan);
h1 = plot(t,Code45(:,4));
h2 = plot(t, yvals2, 'o', 'MarkerEdgeColor', 'magenta', 'LineStyle', 'none');
legend([h1,h2], "Fit", "Auth data");
hold off

%% Difference between Auth and Base
base_diff = transpose(yvals2) - Code45(:,4);
plot(t, base_diff)
title("Author - My prediction")

%% Looping through glycopermeability to see how data changes
glyco_perm = linspace(1,2,50);
tspan = linspace(0,2300,2400);
y0 = [5.73e-3,8e-4,4e-3,4.7e-3]; %femto-mol/um^3

%         Time point   Compartment     Permeability value
M = zeros(length(tspan), 4, length(glyco_perm));
for i = 1:length(glyco_perm)
    theta_loop = [glyco_perm(i), 1.68*10^3];
    [t,Code45] = ode45( @(t,y)matcal_system(t,y,theta_loop), tspan, y0);
    M(:,:,i) = Code45;
end

%% Plotting the results of permeability variation
% Decreasing permeability barely affects how much glucose is in the
% glycocalyx. Slight decrease, but barely anything tbh. 
% Extract the data for plotting
data = squeeze(M(:, 2, :));
% Use the parula colormap
cmap = parula(size(data, 2));
% Number of lines to plot
numLines = size(data, 2);
% Create the plot
figure(1);
hold on;
% Plot each line with a different color from the colormap
for i = 1:numLines
    plot(t, data(:, i), 'Color', cmap(i, :), 'LineWidth', 2);
end
title('Glycocalyx (blue->yellow)')
hold off;

% Decreasing permeability I think barely affects the endothelial cells 
% Extract the data for plotting
data = squeeze(M(:, 3, :));
% Use the parula colormap
cmap = parula(size(data, 2));
% Number of lines to plot
numLines = size(data, 2);
% Create the plot
figure(2);
hold on;
% Plot each line with a different color from the colormap
for i = 1:numLines
    plot(t, data(:, i), 'Color', cmap(i, :), 'LineWidth', 2);
end
ylim([0,3e-5])
title('EC (blue->yellow)')
hold off;

% Adipose tissue however does change quite a bit 
data = squeeze(M(:, 4, :));
% Use the parula colormap
cmap = parula(size(data, 2));
% Number of lines to plot
numLines = size(data, 2);
% Create the plot
figure(3);
hold on;
% Plot each line with a different color from the colormap
for i = 1:numLines
    plot(t, data(:, i), 'Color', cmap(i, :), 'LineWidth', 2);
end
title('Adipose (blue->yellow)')
hold off;

%% Glucose in adipose tissue at 500s
M(500,4,:)


%% Looping through glyco volume to see how data changes
glyco_thickness = linspace(0.001,0.1,50);
tspan = linspace(0,300,400);
y0 = [5.73e-3,8e-4,4e-3,4.7e-3]; %femto-mol/um^3

%         Time point   Compartment     Permeability value
M = zeros(length(tspan), 4, length(glyco_thickness));
for i = 1:length(glyco_thickness)
    glyco_vol = 600 * pi * ((5-0.5)^2 - (5-0.5-glyco_thickness(i))^2);
    theta_loop = [1.4064, glyco_vol];
    [t,Code45] = ode45( @(t,y)matcal_system(t,y,theta_loop), tspan, y0);
    M(:,:,i) = Code45;
end

%% Plotting the results of volume variation
% Decreasing permeability barely affects how much glucose is in the
% glycocalyx. Slight decrease, but barely anything tbh. 
% Extract the data for plotting
data = squeeze(M(:, 2, :));
% Use the parula colormap
cmap = parula(size(data, 2));
% Number of lines to plot
numLines = size(data, 2);
% Create the plot
figure(1);
hold on;
% Plot each line with a different color from the colormap
for i = 1:numLines
    plot(t, data(:, i), 'Color', cmap(i, :), 'LineWidth', 2);
end
title('Glycocalyx (blue->yellow)')
hold off;

% Decreasing permeability I think barely affects the endothelial cells 
% Extract the data for plotting
data = squeeze(M(:, 3, :));
% Use the parula colormap
cmap = parula(size(data, 2));
% Number of lines to plot
numLines = size(data, 2);
% Create the plot
figure(2);
hold on;
% Plot each line with a different color from the colormap
for i = 1:numLines
    plot(t, data(:, i), 'Color', cmap(i, :), 'LineWidth', 2);
end
ylim([0,3e-5])
title('EC (blue->yellow)')
hold off;

% Adipose tissue however does change quite a bit 
data = squeeze(M(:, 4, :));
% Use the parula colormap
cmap = parula(size(data, 2));
% Number of lines to plot
numLines = size(data, 2);
% Create the plot
figure(3);
hold on;
% Plot each line with a different color from the colormap
for i = 1:numLines
    plot(t, data(:, i), 'Color', cmap(i, :), 'LineWidth', 2);
end
title('Adipose (blue->yellow)')
hold off;

%% Functions
function Cdot = matcal_system(~,y,theta)
    % [vr, vg, ve, va, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10];
    vr = 1.06*10^4; % um^3 volume
    vg = theta(2);
    ve= 8.95*10^3;
    va= 4.82*10^5;

    k1=0;
    k2=0;
    % k3= 1.4064; %fmol/s % increasing this I believe is increasing permeability
    k4=149.9983;
    k5=0.3438; % 
    k6=0;
    k7=35.0075; % 
    k8=1.11;
    k9= 5.9959e-3;% 1.032*5.9959e-3 * 0.7;% 5.78/(va*2e-3); % nconsumption
    k10=0;

    eq1 = 0; %(1/vr)*(-k1*vr*y(1) - k2*(vr*y(1))/(va*y(4)) - k3*vr*y(1));  
    eq2 = (1/vg)*( theta(1)*vr*y(1) - k4*vg*y(2) - k5*vg*y(2) + k6*ve*y(3) - k7*vg*y(2) + k10*va*y(4) );
    eq3 = (1/ve)*( k5*vg*y(2) - k6*ve*y(3) - k8*ve*y(3) );
    eq4 = (1/va)*(-k9*va*y(4) + k8*ve*y(3) + k7*vg*y(2) - k10*va*y(4) + k2*vr*y(1) );
    %             coefficient*-5.7837
    Cdot = [eq1;eq2;eq3;eq4];
end