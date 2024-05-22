% MPS2_improvement.m
% Source: https://www.sciencedirect.com/science/article/pii/S0022519321003027
% Purpose:
% Take the author's simple model of transport across the endothelial
% barrier and improve upon it with the details of a glycocalyx layer.
% Authors: Matthew Kozubov, Cal McKinney

% intial conditions and ode solver
% artery, glyco, EC, adipo
C0 = [6e-3;8e-4;4e-3;4.7e-3]; %femto-mol/um^3
tspan = [0 10];% [0 0.01];
opts = odeset('RelTol',1e-5,'AbsTol',1e-7);
tiledlayout(2,3);

% Testing setup to change one variable at a time
coefficient = [0;0.01;1;2;100;1000]; %[0;0.001;0.003;0.006;0.009;0.01]; %
for i = 1:length(coefficient)
    sprintf("Plotting %i", i)

    % Running simulation
    [t,Code45] = ode45( @(t,y)matcal_system(t,y,coefficient(i)), tspan, C0, opts);
    
    % Plotting
    nexttile
    hold on
    plot(t,Code45)
    ylabel("Concentration of glucose (fmol/um^3)")
    xlabel("Time (seconds)")
    ylim([0,6.5e-3])
    title(sprintf("\nGlucose concentrations with coef %.3f\n", coefficient(i)))
    legend(["Artery", "Glyco", "EC", "Tissue"])
end


function Cdot = matcal_system(t,y,coefficient)
    % Transport and volume constants
    vr = 1.06*10^4; % um^3 volume
    vg = 1.68*10^3;
    ve= 8.95*10^3;
    va= 4.82*10^5;
    
    k1=0;
    k2=0;
    k3=coefficient*0.0958; % (6.105)/(2*vr*6e-3); %fmol/s % increasing this I believe is increasing permeability
    k4=0;
    k5=0.3363; % 0.452 / (vg*8e-4)
    k6=0;
    k7=35;%(6.105*va*2e-3)/(2*vg*8e-4);
    k8=1.11;
    k9=5.9959e-3;% 5.78/(va*2e-3); % nconsumption
    k10=0;

    eq1 = 0; %(1/vr)*(-k1*vr*y(1) - k2*(vr*y(1))/(va*y(4)) - k3*vr*y(1));  
    eq2 = (1/vg)*( k3*vr*y(1) - k4*vg*y(2) - k5*vg*y(2) + k6*ve*y(3) - k7*(vg*y(2))/(va*y(4)) + k10*(va*y(4))/(vg*y(2)) );
    eq3 = (1/ve)*( k5*vg*y(2) - k6*ve*y(3) - k8*(ve*y(3))/(va*y(4)) );
    eq4 = (1/va)*( -5.7837 + k8*(ve*y(3))/(va*y(4)) + k7*((vg*y(2))/(va*y(4))) - k10*(va*y(4))/(vg*y(2)) + k2*(vr*y(1))/(va*y(4)) );
    %              -k9*va*y(4)
    Cdot = [eq1;eq2;eq3;eq4];
end






