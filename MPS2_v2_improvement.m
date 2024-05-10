% MPS2_improvement.m
% Source: https://www.sciencedirect.com/science/article/pii/S0022519321003027
% Purpose:
% Take the author's simple model of transport across the endothelial
% barrier and improve upon it with the details of a glycocalyx layer.
% Authors: Matthew Kozubov, Cal McKinney

% intial conditions and ode solver
C0 = [100;50;5;10;20]; % [Cr, Cv, Cg, Ce, Ca]
tspan = [0 5]
opts = odeset('RelTol',1e-5,'AbsTol',1e-7);
[t,Code45] = ode45(@matcal_system, tspan, C0, opts);

plot(t,Code45)
ylabel("Concentration of glucose (mol/L)")
xlabel("Time (seconds)")
title(sprintf("\nPreliminary plot of glucose concentrations based on made up data\n"))
legend(["Artery", "Venous", "Glyco", "EC", "Tissue"])


function Cdot = matcal_system(t,y)
    % Transport and volume constants
    vr=1;
    vv=1;
    vg=0.1;
    ve=0.3;
    va=1;
    vdot1=0.5;
    vdot2=0.05;
    vdot3=0.1;
    vdot4=0.05;
    nk5=20;
    nk6=1;
    nk7=20;
    nk8=3;
    nk9=20;
    nk10=0.2;
    constants = [vr;vv;vg;ve;va;
                vdot1;vdot2;vdot3;vdot4;
                nk5;nk6;
                nk7;nk8;
                nk9;nk10];

    eq1 = (-1/vr)*((vdot1+vdot2+vdot3)*y(1));
    eq2 = (1/vv)*(vdot1*y(1) + vdot4*y(3));
    eq3 = (1/vg)*(vdot3*y(1) - vdot4*y(3) + nk6 + nk10 - nk5 - nk7);
    eq4 = (1/ve)*(nk5 - nk6 - nk8);
    eq5 = (1/va)*(vdot2*y(1) + nk7 + nk8 - nk9 - nk10);
    Cdot = [eq1;eq2;eq3;eq4;eq5];
end






