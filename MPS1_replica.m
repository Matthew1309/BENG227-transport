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

% % Definition of variables
% mdot1 = ;
% mdot2 = ;
% mdot3 = 9.52e-5; % m^3/s
% mdot4 = ;
% 
% CG1 = ;
% CG2 = ;
% CG3 = ;
% 
% rho4 = rho3 = rhob = 1050; % kg/m^3
% 
% Vdot1 = ;
% Vdot2 = ;
% Vdot3 = ;
% g = 9.8; % m/s^2
% tcp = 7.32; % s
% Cap_vessel_diam = ;
% 
% artery_size_M = readmatrix("phen model params - Sheet1.csv", Range=[2, 4])
% 
% %%%%%%%%
% %%% PS1
% %%%%%%%%
% % Flow out of right atrium is equal to inflow with vena inflows
% mdot3 = mdot1 + mdot2; 
% CG3 = (CG1.*Vdot1 + CG2.*Vdot2) ./ Vdot3;
% 
% %%%%%%%%
% %%% PS2 and 3
% %%%%%%%%
% mdot4 = mdot3;
% CG4 = CG3;
% 
% % Mechanical energy balance
% P3./rho3 + g.*z3 + (v3.^2)/2 + nabla .* What = P4 ./ rho4 + g.*z4 + (v4.^2)/2;
% What = (P4 - P3) ./ rhob;
% 
% %%%%%%%%
% %%% PS4
% %%%%%%%%
% mdot7 = mdot4;
% mdot6 = mdot5;
% 
% CG5 = CG6 = 0;
% Vdot4 = Vdot7; 
% CG7 = CG4;
% 
% %%%%%%%%
% %%% PS5 and 6
% %%%%%%%%
% mdot8 = mdot7; 
% CG8 = CG7; 
% What = (P8 - P7)./rhob;
% 
% %%%%%%%%
% %%% PS7 and 8
% %%%%%%%%
% % Pretend only one vessel
% Vdot9input = Vdot9output + Vdot10output;
% Vdot9output = PC .* Vdot9input;
% Vdot10output = (1-PC) .* Vdot9input;
% 
% CG9 = CG9 = CG10 ; %
% % Big vessels
% P9out = P9in + (rhob.*(v9in.^2 - v9out.^2))./2 - rhob .* h .* finout;
% % Small vessels
% P9out = P9in - (128.*mu9.*L9.*Vdot9)./(pi.*Cap_vessel_diam.^4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pressure monstrousity

%%%%%%
% Constants
%%%%%%
%% Initialize constants
format longg
artery_size_M = readmatrix("phen model params - Sheet1.csv", Range=[2, 4]);
D_9_in_M = artery_size_M(:,1) .* 10^-3; % meter units
D_9_out_M = artery_size_M(:,2) .* 10^-3;
D_10_out_M = artery_size_M(:,3) .* 10^-3;
vessel_length_M = artery_size_M(:,4) .* 10^-3;

blood_density = 1060; %kg/m^3;
H = 0.45;
epsilon = 1e4;
K1 = [500; 150; 1000]; % all of them are y-valves
Kinf = [0.15; 0.5; 2]; % all of them are y-valves

%%%%%%
% Testing
%%%%%%
%% Climb the xs area branch
% Compute xs areas
CSA_9_in_M = cross_sectional_area(D_9_in_M);
CSA_9_out_M = cross_sectional_area(D_9_out_M);
CSA_10_out_M = cross_sectional_area(D_10_out_M);

% Compute partition coefficients
PC9 = partion_coeff_9(CSA_9_out_M, CSA_10_out_M);

% Compute volumetric flows
Vol0_9 = 9.52*10^-5; % m8/blood_density (but not really, just m3) I think they messed up 
Vol_9 = zeros(length(PC9),1);
last_flow = Vol0_9;
for i = 1:(length(PC9))
    % fprintf('%i\n', i)
    Vol_9(i) = volumetric_flow(PC9(i), last_flow);
    last_flow = Vol_9(i);
end

% Compute velocity in tube
velocity_9in = velocity(Vol_9, CSA_9_in_M); % m/s units
velocity_10out = velocity(Vol_9, CSA_10_out_M);
velocity_9out = velocity(Vol_9, CSA_9_out_M);

%% Climb the Bj branch
Bjs = Bj(D_9_in_M .* 10^3); % Expects input in units of mm;
Sjs = Sj(Bjs, D_9_in_M .* 10^3);
apparent_vessel_visc_j = apparent_vessel_viscosity(D_9_in_M .* 10^3);
vessel_viscosity_j = vessel_viscosity(apparent_vessel_visc_j, Sjs, H);
reynolds_number_j = reynolds_number(vessel_viscosity_j, velocity_9in, D_9_in_M .* 10^3, blood_density);

%% Climb the D'', Darcy, Kexp/Kcon branches
dapostrophie_j = dapostrophie(D_9_in_M.*10^3);
Kacc_j = zeros(length(dapostrophie_j),3);
for i = 1:length(dapostrophie_j)
    Kacc_j(i,:) = Kacc(reynolds_number_j(i), dapostrophie_j(i), K1, Kinf);
end

darcy_friction_coeff_j = darcy_friction_coeff(reynolds_number_j, D_9_in_M .* 10^3, epsilon);

Kexp_j = Kexp(D_9_out_M .* 10^3, D_10_out_M .* 10^3);
Kcon_j = Kcon(D_9_out_M .* 10^3, D_10_out_M .* 10^3);

friction_loss_j = friction_loss(darcy_friction_coeff_j, vessel_length_M .* 10^3, velocity_9out, velocity_10out, Kacc_j, Kexp_j, Kcon_j)

%% Pressure computation
pressure_out_j = zeros(length(PC9),1);
P0_9 = 94.23; %%mmHg
last_pressure = P0_9;
for i = 1:(length(PC9))
    pressure_out_j(i) = pressure_out_large_vessel(last_pressure, velocity_9in(i), velocity_9out(i), friction_loss_j(i), blood_density);
    last_flow = pressure_out_j(i);
end
pressure_out_j

%%%%%%
% Functions 
%%%%%%
%%
function output = pressure_out_large_vessel(pressure_in, velocity9_in, velocity9_out, friction_loss, blood_density)
    % Representing one vessel, pass this
    % pressure on the input 
    % velocity on the input and output that matter
    % loss of pressure due to friction
    % it will output one value representing the pressure on the output
    output = pressure_in + (blood_density.*(velocity9_in.^2 - velocity9_out.^2)) ./ 2 - (blood_density .* friction_loss);
end

function velocity = velocity(vol_flow_rate, cross_area)
    % Representing the velocity in one vessel
    % given the volumetric flow rate in this vessel
    % and its cross-sectional area, we can compute the 
    % velocity of the liquid. General function for ins and outs
    velocity = vol_flow_rate ./ cross_area;
end

function friction_loss = friction_loss(darcy_friction_coeff, vessel_length, velocity9_out, velocity10_out, Kacc, Kexp, Kcon)
    % Passing the darcy coeff for the vessel, the length, velocities, and
    % several coefficients we can compute the friction loss. Additionally
    % pass a 3-vector of Kacc terms for the vessel computed in a different
    % function.
    middle_term = zeros(length(vessel_length),1);
    for i = 1:length(vessel_length)
        middle_term(i) = sum( Kacc(i,:) .* max([velocity9_out(i), velocity10_out(i)]) );
    end
    friction_loss = (darcy_friction_coeff .* vessel_length .* velocity9_out) + (middle_term./2) + (Kexp .* velocity9_out)./2 + (Kcon.*velocity9_out)./2;
end

function darcy_friction_coeff = darcy_friction_coeff(reynolds_number, diameter9_in, epsilon)
    if reynolds_number <= 2000
        darcy_friction_coeff = 64./reynolds_number;
    else
        term1 = epsilon ./ (diameter9_in./3.71);
        coef2 = 502./reynolds_number;
        term3 = 14.5./reynolds_number;
        darcy_friction_coeff = (-2.*log(term1 - coef2.*(term1 + term3)) ).^-2;
    end
end

function Kacc = Kacc(reynolds_number, dapostrophie, K1, Kinf)
    Kacc = zeros(3,1);
    for i = 1:3
        Kacc(i) = K1(i)./reynolds_number + Kinf(i).*(1+(1./dapostrophie));
    end
end

function dapostrophie = dapostrophie(vessel_diameter)
    dapostrophie = 39.37 .* vessel_diameter;
end

function Kexp = Kexp(vessel_diameter_in, vessel_diameter_out)
    Kexp = (1-(vessel_diameter_in ./ vessel_diameter_out).^2).^2;
end

function Kcon = Kcon(vessel_diameter_in, vessel_diameter_out)
    Kcon = 0.5.*(1-(vessel_diameter_out ./ vessel_diameter_in).^2).^2;
end

function area = cross_sectional_area(vessel_diameter)
    area = 0.25 .* (pi .* vessel_diameter.^2);
end

function reynolds_number = reynolds_number(vessel_viscosity, velocity_9_in, vessel_diameter_in, blood_density)
    reynolds_number = (blood_density.*velocity_9_in.*vessel_diameter_in) ./ vessel_viscosity;
end

function vessel_viscosity = vessel_viscosity(apparent_vessel_visc, Sj, H)
    vessel_viscosity = 1 + ((apparent_vessel_visc - 1) .* (((1-H).^Sj)./(0.55.^Sj)) );
end

function apparent_vessel_viscosity = apparent_vessel_viscosity(vessel_diameter)
    term1 = 220 .* exp(-1.3 .* vessel_diameter .* 10^-6);
    term2 = 2.44 .* exp(-0.06 .* (vessel_diameter .* 10^-6).^0.645);
    apparent_vessel_viscosity = 3.2 + term1 - term2;
end

function S = Sj(Bj, vessel_diameter_in)
    term1 = 0.8 + exp(-0.075 .* vessel_diameter_in .* 10^-6);
    S = (term1 .* (Bj-1)) + Bj;
end

function B = Bj(vessel_diameter_in)
    % Original paper wants diameter in units of micrometers
    % Our authors I believe already converted to units of mm implicityly
    B = (1 + (vessel_diameter_in.^12 .* 10^-6)).^-1;
end

function vol_flow = volumetric_flow(partion_coeff, flow_coming_in)
    vol_flow = partion_coeff * flow_coming_in;
end

function PC = partion_coeff_9(cross_sec_9_out, cross_sec_10_out)
    PC = cross_sec_9_out ./ (cross_sec_9_out + cross_sec_10_out);
end
