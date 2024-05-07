% MPS2_improvement.m
% Source: https://www.sciencedirect.com/science/article/pii/S0022519321003027
% Purpose:
% Take the author's simple model of transport across the endothelial
% barrier and improve upon it with the details of a glycocalyx layer.
% Authors: Matthew Kozubov, Cal McKinney

% intial conditions and ode solver
C0 = [100;0;0;0];
tspan = [0 5]
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
    vdotr=0.5;
    vdotg=0.5;
    nintocell=1;
    noutofcellup=1;
    nthroughjunc=30;
    noutofcelldown=3;
    nconsumption=20;
    constants = [vr;vg;ve;va;
                vdotr;vdotg;
                nintocell;noutofcellup;
                nthroughjunc;noutofcelldown;
                nconsumption];

    eq1 = (-1/vr)*(vdotr*y(1));
    eq2 = (1/vg)*(vdotr*y(1) - vdotg*(y(2)) - nintocell + noutofcellup + nthroughjunc);
    eq3 = (1/ve)*(nintocell - noutofcellup - noutofcelldown);
    eq4 = (1/va)*(noutofcelldown + nthroughjunc - nconsumption);
    Cdot = [eq1;eq2;eq3;eq4];
end






