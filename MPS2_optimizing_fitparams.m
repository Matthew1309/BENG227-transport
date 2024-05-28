% MPS2_improvement.m
% Source: https://www.sciencedirect.com/science/article/pii/S0022519321003027
% Purpose:
% Take the author's simple model of transport across the endothelial
% barrier and improve upon it with the details of a glycocalyx layer.
% Authors: Matthew Kozubov, Cal McKinney
% This version almost has it mimicking the author's data, but not entirely

author_fig16 = readmatrix("phenomenological-based-model-adipose-tissue-glucose-levels.csv");
author_fig16 = author_fig16(2:end,:);
author_fig16(:,1) = author_fig16(:,1) - author_fig16(1,1);
plot(author_fig16(:,1), author_fig16(:,2));


% Transport and volume constants
k1=0;
k2=0;
k3=2*0.0958; % (6.105)/(2*vr*6e-3); %fmol/s % increasing this I believe is increasing permeability
k4=0;
k5=0.3363; % 0.452 / (vg*8e-4)
k6=0;
k7=35;%(6.105*va*2e-3)/(2*vg*8e-4);
k8=1.11;
k9=1.032*5.9959e-3 * 0.7;% 5.78/(va*2e-3); % nconsumption
k10=0;
% theta = [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10];
theta = [k3, k7];


% intial conditions and ode solver
% artery, glyco, EC, adipo
y0 = [5.73e-3,8e-4,4e-3,4.7e-3]; %femto-mol/um^3
tspan = linspace(0,2450);
opts = odeset('RelTol',1e-5,'AbsTol',1e-7);
% How to pass theta to ode45
[t,Code45] = ode45( @(t,y)matcal_system(t,y,theta), tspan, y0, opts);
yvalstrue = interp1(author_fig16(:,1),author_fig16(:,2),tspan);
yvalstrue = [repelem(6.3e-3,length(yvalstrue));...
             repelem(0.1e-3,length(yvalstrue));...
             repelem(2e-3,length(yvalstrue));...
             yvalstrue.*10^-3];

%% 
% Following tutorial for optimizing params if I know all my variables
% https://www.mathworks.com/help/optim/ug/fit-ode-problem-based-least-squares.html
r = optimvar('theta', 2, 'LowerBound', 1e-2, 'UpperBound', 300);
myfcn = fcn2optimexpr(@RtoODE,r,tspan,y0);
obj = sum(sum(myfcn - yvalstrue).^2);
prob = optimproblem("Objective",obj);
r0.theta = theta;
[rsol,sumsq] = solve(prob,r0)
save('myData.mat', 'rsol', 'sumsq');

[t,Code45] = ode45( @(t,y)matcal_system(t,y,rsol.theta), tspan, y0, opts);
plot(t,Code45)

%% 
% % Now I only know the adipose tissue and am optimizing for that only
r = optimvar('theta', 2, 'LowerBound', 1e-2, 'UpperBound', 300);
myfcn2 = fcn2optimexpr(@RtoODE2, r, tspan, y0);
yvals2 = yvalstrue([4],:);
obj2 = sum(sum(myfcn2 - yvals2).^2);
prob2 = optimproblem("Objective",obj2);
r0.theta = theta;
[rsol2,sumsq2] = solve(prob2,r0)
rsol2.theta

%%
[t,Code45] = ode45( @(t,y)matcal_system(t,y,rsol2.theta), tspan, y0, opts);
figure(1)
hold on
h1 = plot(t,Code45(:,1));
h3 = plot(t,Code45(:,2));
h4 = plot(t,Code45(:,3));
h5 = plot(t,Code45(:,4));
h2 = plot(t, yvals2, 'o', 'MarkerEdgeColor', 'magenta', 'LineStyle', 'none');
legend([h1,h3,h4,h5,h2], ["Artery", "Glyco", "EC", "Adipose", "Auth data"]);
hold off

figure(2)
hold on
h1 = plot(t,Code45(:,4));
h2 = plot(t, yvals2, 'o', 'MarkerEdgeColor', 'magenta', 'LineStyle', 'none');
legend([h1,h2], "Fit", "Auth data");
hold off

%% Functions
function Cdot = matcal_system(~,y,theta)
    % [vr, vg, ve, va, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10];
    vr = 1.06*10^4; % um^3 volume
    vg = 1.68*10^3;
    ve= 8.95*10^3;
    va= 4.82*10^5;

    k1=0;
    k2=0;
    % k3=2*0.0958; % (6.105)/(2*vr*6e-3); %fmol/s % increasing this I believe is increasing permeability
    k4=0;
    k5=0.3363; % 0.452 / (vg*8e-4)
    k6=0;
    % k7=35;%(6.105*va*2e-3)/(2*vg*8e-4);
    k8=1.11;
    k9=1.032*5.9959e-3 * 0.7;% 5.78/(va*2e-3); % nconsumption
    k10=0;

    eq1 = 0; %(1/vr)*(-k1*vr*y(1) - k2*(vr*y(1))/(va*y(4)) - k3*vr*y(1));  
    eq2 = (1/vg)*( theta(1)*vr*y(1) - k4*vg*y(2) - k5*vg*y(2) + k6*ve*y(3) - theta(2)*vg*y(2) + k10*va*y(4) );
    eq3 = (1/ve)*( k5*vg*y(2) - k6*ve*y(3) - k8*ve*y(3) );
    eq4 = (1/va)*(-k9*va*y(4) + k8*ve*y(3) + theta(2)*vg*y(2) - k10*va*y(4) + k2*vr*y(1) );
    %             coefficient*-5.7837
    Cdot = [eq1;eq2;eq3;eq4];
end

function solpts = RtoODE(theta,tspan,y0)
    sol = ode45( @(t,y)matcal_system(t,y,theta), tspan, y0);
    solpts = deval(sol,tspan);
end

function solpts = RtoODE2(theta,tspan,y0)
    sol = ode45( @(t,y)matcal_system(t,y,theta), tspan, y0);
    solpts = deval(sol,tspan);
    solpts = solpts(4,:); % Just the adipose tissue
end