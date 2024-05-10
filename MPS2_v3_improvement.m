% MPS2_improvement.m
% Source: https://www.sciencedirect.com/science/article/pii/S0022519321003027
% Purpose:
% Take the author's simple model of transport across the endothelial
% barrier and improve upon it with the details of a glycocalyx layer.
% Authors: Matthew Kozubov, Cal McKinney

% intial conditions and ode solver
C0 = [100;2;2;20];
tspan = [0 10];
opts = odeset('RelTol',1e-5,'AbsTol',1e-7);
[t,Code45] = ode45(@matcal_system, tspan, C0, opts);

plot(t,Code45)
ylabel("Concentration of glucose (mol/L)")
xlabel("Time (seconds)")
title(sprintf("\nPreliminary plot of glucose concentrations based on made up data\n"))
legend(["Artery", "Glyco", "EC", "Tissue"])


function Cdot = matcal_system(t,y)
    % Transport and volume constants
    vr = 1;
    vg = 0.1;
    ve=0.3;
    va=1;
    
    k8 = 0.1;
    vdot1 = 0.25;
    k2 = 0.25;
    vdot3 = 0.2;
    vdot4 = 0.15;
    vdot5 = 0.5;
    vdot6 = 0.01;
    k7 = 3*k8;
    vdot7 = k7*vg;
    nconsumption = 5;
    k10 = k8*0.5;

    eq1 = (1/vr)*(-vdot1*y(1) - k2*(vr*y(1))/(va*y(4)) - vdot3*y(1));
    eq2 = (1/vg)*( vdot3*y(1) - vdot4*y(2) - vdot5*y(2) + vdot6*y(3) - vdot7*y(2) + k10*(va*y(4))/(vg*y(2)) );                         ; 
    eq3 = (1/ve)*( vdot5*y(2) - vdot6*y(3) - k8*(ve*y(3))/(va*y(4)) );
    eq4 = (1/va)*( -nconsumption + k8*(ve*y(3))/(va*y(4)) + vdot7*y(2) - k10*(va*y(4))/(vg*y(2)) + k2*(vr*y(1))/(va*y(4)) );
    Cdot = [eq1;eq2;eq3;eq4];
end






