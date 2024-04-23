% MPS1_replica.m
% Source: https://www.sciencedirect.com/science/article/pii/S0022519321003027
% Purpose:
% Replicate the model of glucose circulation proposed by this group
% and include a new level of compartmentalization for diffusion of
% glucose out of the capillaries. The authors didn't include the glycocalyx
% opting instead to model the endothelial layer as just a monolith. Cal and I
% can split endothelial cells into a glycocalyx layer and an endothelial
% layer to show how glucose actually travels across. Might not change
% modeling efficiency, but it can provide insights into the inner workings
% of the system.
% Authors: Matthew Kozubov, Cal McKinney

% Definition of variables
mdot1 = ;
mdot2 = ;
mdot3 = 9.52e-5; % m^3/s
mdot4 = ;

CG1 = ;
CG2 = ;
CG3 = ;

rho4 = rho3 = rhob = 1050; % kg/m^3

Vdot1 = ;
Vdot2 = ;
Vdot3 = ;
g = 9.8; % m/s^2
tcp = 7.32; % s
Cap_vessel_diam = ;

%%%%%%%%
%%% PS1
%%%%%%%%
% Flow out of right atrium is equal to inflow with vena inflows
mdot3 = mdot1 + mdot2; 
CG3 = (CG1.*Vdot1 + CG2.*Vdot2) ./ Vdot3;

%%%%%%%%
%%% PS2 and 3
%%%%%%%%
mdot4 = mdot3;
CG4 = CG3;

% Mechanical energy balance
P3./rho3 + g.*z3 + (v3.^2)/2 + nabla .* What = P4 ./ rho4 + g.*z4 + (v4.^2)/2;
What = (P4 - P3) ./ rhob;

%%%%%%%%
%%% PS4
%%%%%%%%
mdot7 = mdot4;
mdot6 = mdot5;

CG5 = CG6 = 0;
Vdot4 = Vdot7; 
CG7 = CG4;

%%%%%%%%
%%% PS5 and 6
%%%%%%%%
mdot8 = mdot7; 
CG8 = CG7; 
What = (P8 - P7)./rhob;

%%%%%%%%
%%% PS7 and 8
%%%%%%%%
% Pretend only one vessel
Vdot9input = Vdot9output + Vdot10output;
Vdot9output = PC .* Vdot9input;
Vdot10output = (1-PC) .* Vdot9input;

CG9 = CG9 = CG10 ; %
% Big vessels
P9out = P9in + (rhob.*(v9in.^2 - v9out.^2))./2 - rhob .* h .* finout;
% Small vessels
P9out = P9in - (128.*mu9.*L9.*Vdot9)./(pi.*Cap_vessel_diam.^4);

