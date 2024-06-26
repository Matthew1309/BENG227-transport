% MPS2_improvement.m
% Source: https://www.sciencedirect.com/science/article/pii/S0022519321003027
% Purpose:
% Take the author's simple model of transport across the endothelial
% barrier and improve upon it with the details of a glycocalyx layer.
% Authors: Matthew Kozubov, Cal McKinney
% This version almost has it mimicking the author's data, but not entirely

%% Read in the author data
global author_fig16;
author_fig16 = readmatrix("phenomenological-based-model-adipose-tissue-glucose-levels.csv");
author_fig16 = author_fig16(2:end,:);
author_fig16(:,1) = author_fig16(:,1) - author_fig16(1,1);
plot(author_fig16(:,1), author_fig16(:,2));

%% Making a function that can use our optimized params
close all
%       1.4064
theta = [1.4064, 1.68*10^3];

figure(1)
fig = gcf;
fig.Position = [100, 100, 800, 800]; % [left, bottom, width, height] in pixels
tiledlayout(2,1);

nexttile
hold on
y0 = [5.73e-3,8e-4,4e-3,4.7e-3]; %femto-mol/um^3
tspan = linspace(0,2450);
[t,Code45] = ode45( @(t,y)matcal_system(t,y,theta), tspan, y0);
h1 = plot(t,Code45(:,1), 'LineWidth', 2);
h3 = plot(t,Code45(:,2), 'LineWidth', 2);
h4 = plot(t,Code45(:,3), 'LineWidth', 2);
h5 = plot(t,Code45(:,4), 'LineWidth', 2);
h2 = plot(t, author_y_vals(t), 'o', "MarkerSize", 2 ,'MarkerEdgeColor', 'magenta', 'LineStyle', 'none');
% legend([h1,h3,h4,h5,h2], ...
  %   ["Artery", "Glyco", "EC", "Adipose", "Auth data"], ...
    % "Location", 'east');
title(sprintf('Predicted glucose concentrations in compartments\nwith author optimized values (long time scale)'))
% Add labels and title
xlabel('Time (s)');
ylabel('Glucose Conc. (mmol/L)'); 
% Add grid lines
grid off;
% Customize tick marks and labels
set(gca, 'FontSize', 12); % Adjust font size
% Adjust plot appearance
set(gcf, 'Color', 'w'); % Set background color to white
set(gca, 'Box', 'on'); % Draw box around the plot
legend([h1,h3,h4,h5,h2], ...
    ["Artery", "Glyco", "EC", "Adipose", "Auth data"], ...
    'Location', 'eastoutside');
hold off

nexttile
hold on
y0 = [5.73e-3,8e-4,4e-3,4.7e-3]; %femto-mol/um^3
tspan = linspace(0,10);
[t,Code45] = ode45( @(t,y)matcal_system(t,y,theta), tspan, y0);
h1 = plot(t,Code45(:,1), 'LineWidth', 2);
h3 = plot(t,Code45(:,2), 'LineWidth', 2);
h4 = plot(t,Code45(:,3), 'LineWidth', 2);
h5 = plot(t,Code45(:,4), 'LineWidth', 2);
h2 = plot(t, author_y_vals(t), 'o', "MarkerSize", 2 ,'MarkerEdgeColor', 'magenta', 'LineStyle', 'none');
title(sprintf('Predicted glucose concentrations in compartments\nwith author optimized values (short time scale)'))
% Add labels and title
xlabel('Time (s)');
ylabel('Glucose Conc. (mmol/L)'); 
% Add grid lines
grid off;
% Customize tick marks and labels
set(gca, 'FontSize', 12); % Adjust font size
% Adjust plot appearance
set(gcf, 'Color', 'w'); % Set background color to white
set(gca, 'Box', 'on'); % Draw box around the plot
hold off
legend([h1,h3,h4,h5,h2], ...
    ["Artery", "Glyco", "EC", "Adipose", "Auth data"], ...
    'Location', 'eastoutside');
set(gca, 'Position', [0.1, 0.1, 0.7, 0.9]);

% Save the plot as a high-resolution image for publication
print('long_and_short_timescale_all_compartment_concentrations.png','-dpng','-r300'); % Save as PNG with 300 DPI

figure(2)
fig = gcf;
fig.Position = [100, 100, 900, 400]; % [left, bottom, width, height] in pixels

hold on
tspan = linspace(0,2450);
[t,Code45] = ode45( @(t,y)matcal_system(t,y,theta), tspan, y0);
yvals2 = 10^-3 * interp1(author_fig16(:,1),author_fig16(:,2),tspan);
h1 = plot(tspan,Code45(:,4));
h2 = plot(tspan, yvals2, 'o', 'MarkerEdgeColor', 'magenta', 'LineStyle', 'none');
legend([h1,h2], "Fit", "Auth data", 'Location', 'eastoutside');
title(sprintf('Fit vs Author adipose interstitial \nglucose concentration comparison'))
% Add labels and title
xlabel('Time (s)');
ylabel('Glucose Conc. (mmol/L)'); 
% Add grid lines
grid off;
% Customize tick marks and labels
set(gca, 'FontSize', 12); % Adjust font size
% Adjust plot appearance
set(gcf, 'Color', 'w'); % Set background color to white
set(gca, 'Box', 'on'); % Draw box around the plot

hold off
set(gcf, 'Color', 'w');

% Save the plot as a high-resolution image for publication
print('Adipose_tissue_author_vs_recreation_concentrations.png','-dpng','-r300'); % Save as PNG with 300 DPI


%% Difference between Auth and Base
base_diff = transpose(yvals2) - Code45(:,4);
RMSE_adipose = sqrt(mean((base_diff).^2));
average_author_data = mean(Code45(:,4))

format long; % Display more digits
thing = [t, transpose(yvals2),  Code45(:,4), transpose(yvals2) - Code45(:,4),(transpose(yvals2) - Code45(:,4)).^2]

plot(t, base_diff)
% plot(t, base_diff.^2)
title(sprintf("Difference between author model and replica\n RMSE: %.3s", RMSE_adipose))
% Add labels and title
xlabel('Time (s)');
ylabel('Glucose Conc. (mmol/L)'); 
% Add grid lines
grid off;
% Customize tick marks and labels
set(gca, 'FontSize', 12); % Adjust font size
% Adjust plot appearance
set(gcf, 'Color', 'w'); % Set background color to white
set(gca, 'Box', 'on'); % Draw box around the plot

% Save the plot as a high-resolution image for publication
print('difference_between_auth_and_our_fit_models.png','-dpng','-r300'); % Save as PNG with 300 DPI


%% Looping through glycopermeability to see how data changes
glyco_perm = linspace(0.7,1.4064,50);
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
fig = gcf;
fig.Position = [100, 100, 700, 400]; % [left, bottom, width, height] in pixels

hold on;
% Plot each line with a different color from the colormap
for i = 1:numLines
    plot(t, data(:, i), 'Color', cmap(i, :), 'LineWidth', 2);
end
% Add labels and title
xlabel('Time (s)');
ylabel('Glucose Conc. (mmol/L)');
title(sprintf('Predicted interstitial adipose tissue \nglucose concentration with varying permeability'));
% Add grid lines
grid off;
% Customize tick marks and labels
set(gca, 'FontSize', 12); % Adjust font size

% Adjust plot appearance
set(gcf, 'Color', 'w'); % Set background color to white
set(gca, 'Box', 'on'); % Draw box around the plot
% Add colorbar
cb = colorbar;
% Set colormap limits to match your data range
clim([min(glyco_perm), max(glyco_perm)]);
fig = gcf; % Get current figure handle
fig.Position(3) = fig.Position(3) * 1.15; % Increase width by a factor of 1.5 (adjust as needed)

% Save the plot as a high-resolution image for publication
print('Adipose_tissue_concentration.png','-dpng','-r200'); % Save as PNG with 300 DPI

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

function y = author_y_vals(tspan)
    global author_fig16;
    y = 10^-3 * interp1(author_fig16(:,1),author_fig16(:,2),tspan);
end
