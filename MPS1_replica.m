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
artery_size_M(:,1:4) = artery_size_M(:,1:4) .* 10^-3; % if .* 10^-3 units are in meters else mm
D_9_in_M = artery_size_M(:,1);
D_9_out_M = artery_size_M(:,2);
D_10_out_M = artery_size_M(:,3);

D_9_out_M(4) = 0.0121;
D_10_out_M(4) = 0.0275;

vessel_length_M = artery_size_M(:,4);
vessel_accessories = artery_size_M(:,5);

blood_density = 1060; %kg/m^3;
H = 0.45;
epsilon = 1e4;
% 2-K (Hooper) Method
K1 = containers.Map({1,2,3}, [500,150,1000]); % dictionary that the exel sheet defines
Kinf = containers.Map({1,2,3}, [0.15; 0.5; 2]); % dictionary that the exel sheet defines

%%%%%%
% Testing
%%%%%%
%% Climb the xs area branch (i like this section implementation)
% Compute xs areas
CSA_9_in_M = cross_sectional_area(D_9_in_M); %meter^2
CSA_9_out_M = cross_sectional_area(D_9_out_M);
CSA_10_out_M = cross_sectional_area(D_10_out_M);

% Compute partition coefficients
PC9 = partion_coeff_9(CSA_9_out_M, CSA_10_out_M);
% PC9 = ones(length(CSA_9_out_M),1);
% for i = 2:length(CSA_9_out_M)
%     PC9(i) = partion_coeff_9(CSA_9_out_M(i-1), CSA_10_out_M(i-1));
% end

thing = [transpose(0:48), CSA_9_out_M, CSA_10_out_M, PC9];

% Compute volumetric flows
Vol_9 = zeros(length(PC9),1);
Vol_10 = zeros(length(PC9),1);
Vol_9(1) = 9.52*10^-5; % m8/blood_density (but not really, just m3) I think they messed up m^3/s
last_flow = Vol_9(1);
for i = 2:(length(PC9))
    % fprintf('%i\n', i)
    Vol_9(i) = volumetric_flow(PC9(i), last_flow);
    Vol_10(i) = volumetric_flow(1-PC9(i), last_flow);
    last_flow = Vol_9(i);
end

velocity_9in = velocity(Vol_9, CSA_9_in_M); % m/s
velocity_10out = velocity(Vol_10, CSA_10_out_M);
velocity_9out = zeros(length(CSA_9_out_M), 1);
for i = 1:(length(CSA_9_out_M)-1)
    velocity_9out(i) = velocity(Vol_9(i+1), CSA_9_out_M(i));
end

head([velocity_9in, velocity_9out]) 
thing = head([transpose(0:48), velocity_9in, velocity_9out, Vol_9, PC9, CSA_9_in_M]);
plot(velocity_9in)

time_through_vessels = sum(vessel_length_M ./ ((velocity_9in + velocity_9out) ./2));

%% Climb the Bj branch
Bjs = Bj(D_9_in_M); % Expects input in units of mm
Sjs = Sj(Bjs, D_9_in_M);
apparent_vessel_visc_j = apparent_vessel_viscosity(D_9_in_M .* 10^3);
vessel_viscosity_j = vessel_viscosity(apparent_vessel_visc_j, Sjs, H);

reynolds_number_j = reynolds_number(vessel_viscosity_j, ...
                                    velocity_9in, ...
                                    D_9_in_M.*10^3, ...
                                    blood_density);
head([vessel_viscosity_j, reynolds_number_j]) 

%% Climb the D'', Darcy, Kexp/Kcon branches
dapostrophie_j = dapostrophie(D_9_in_M * 10^3); % mm to inches?
Kacc_j = zeros(length(dapostrophie_j),1);
for i = 1:length(dapostrophie_j)
    Kacc_j(i) = Kacc(reynolds_number_j(i), dapostrophie_j(i), K1, Kinf, vessel_accessories(i));
end

% darcy_friction_coeff_j = darcy_friction_coeff(reynolds_number_j, D_9_in_M, epsilon);
darcy_friction_coeff_j = zeros(length(reynolds_number_j), 1);
for i = 1:length(reynolds_number_j)
    darcy_friction_coeff_j(i) = darcy_friction_coeff(reynolds_number_j(i), D_9_in_M(i), epsilon);
end

%Kexp_j = Kexp(D_9_out_M, D_10_out_M);
%Kcon_j = Kcon(D_9_out_M, D_10_out_M);
Kcon_j = zeros(length(D_9_out_M),1);
Kexp_j = zeros(length(D_9_out_M),1);
for i = 2:length(D_9_out_M)
    Kcon_j(i) = Kcon(D_9_in_M(i-1), D_9_out_M(i-1));
    Kexp_j(i) = Kexp(D_9_in_M(i-1), D_9_out_M(i-1));
end

% vessel_length_M .*10^-3
friction_loss_j = friction_loss(darcy_friction_coeff_j, vessel_length_M, velocity_9out, velocity_10out, Kacc_j, Kexp_j, Kcon_j);
head(friction_loss_j)

%% Pressure computation
pressure_out_j = zeros(length(PC9),1);
P0_9 = 133.32 * 94.23; % Pa % 94.23; %%mmHg
last_pressure = P0_9;
for i = 2:(length(PC9))
    pressure_out_j(i-1) = pressure_out_large_vessel(last_pressure, velocity_9in(i-1), velocity_9out(i-1), friction_loss_j(i-1), blood_density);
    last_flow = pressure_out_j(i-1);
end
pressure_out_j = pressure_out_j ./ 133.32


%% Plotting length of vessel traversed vs velocity
floored_vessel_length = ceil(vessel_length_M*1000);
nx = sum(floored_vessel_length);
blood_velocities = zeros(nx,1);
nstart = 1;
for i = 1:length(floored_vessel_length)
    nspots = floored_vessel_length(i);
    blood_velocities(nstart:nstart + nspots - 1) = velocity_9in(i);
    nstart = nstart + nspots;
end
figure(1)
plot(1:nx, blood_velocities)


floored_vessel_length = ceil(vessel_length_M*1000);
nx = sum(floored_vessel_length);
blood_vol_flow = zeros(nx,1);
nstart = 1;
for i = 1:length(floored_vessel_length)
    nspots = floored_vessel_length(i);
    blood_vol_flow(nstart:nstart + nspots - 1) = Vol_9(i);
    nstart = nstart + nspots;
end

figure(2)
plot(1:nx, blood_vol_flow)


floored_vessel_length = ceil(vessel_length_M*1000);
nx = sum(floored_vessel_length);
blood_pressures = zeros(nx,1);
nstart = 1;
for i = 1:length(floored_vessel_length)
    nspots = floored_vessel_length(i);
    blood_pressures(nstart:nstart + nspots - 1) = pressure_out_j(i);
    nstart = nstart + nspots;
end

figure(3)
plot(1:nx, blood_pressures, ".", 'MarkerSize', 8, "LineStyle", "-")
ylim([80, 100])






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
    % Units on friction loss work out so that this whole things should be
    % in pascal = N/m^2
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
    % vessel_length - m
    % velocity - m/s
    % darcy coeff - unitless or 1/m <- probably has to be for units to makes sense

    middle_term = zeros(length(vessel_length),1);
    for i = 1:length(vessel_length)
        middle_term(i) = Kacc(i) .* max([velocity9_out(i).^2, velocity10_out(i).^2]);
    end
    middle_term(isinf(middle_term)) = 0;
    Kexp(isinf(Kexp)) = 0;
    Kcon(isinf(Kcon)) = 0;
    friction_loss = (darcy_friction_coeff .* vessel_length .* velocity9_out.^2) + ...
                    (middle_term./2) + ...
                    (Kexp .* velocity9_out)./2 + ...
                    (Kcon.*velocity9_out)./2;
    % Term 1 - 1/m * m * m^2/s^2
    % Term 2 - 
    % Should have units of - m^2/s^2
    % friction_loss = [(darcy_friction_coeff .* vessel_length .* velocity9_out.^2), ...
      %               (middle_term./2), ...
        %             (Kexp .* velocity9_out)./2 ,...
          %           (Kcon.*velocity9_out)./2];
     % kg/m^3 * m^2/s^2 = Pascal = N/m^2 = kg/m*s^2
end

function darcy_friction_coeff = darcy_friction_coeff(reynolds_number, diameter9_in, epsilon)
    % I think the units are none?
    if reynolds_number <= 2000
        darcy_friction_coeff = 64./reynolds_number;
    else
        term1 = epsilon ./ (diameter9_in./3.71);
        coef2 = 502./reynolds_number;
        term3 = 14.5./reynolds_number;
        darcy_friction_coeff = (-2.*log(term1 - coef2.*(term1 + term3)) ).^-2;
    end
end

function Kacc = Kacc(reynolds_number, dapostrophie, K1, Kinf, vessel_accessory)
    % Passed a vessel one at a time
    % https://neutrium.net/fluid-flow/pressure-loss-from-fittings-2k-method/
    % 2-K (Hooper) Method
    % dapostrophie should be in inches
    Kacc = K1(vessel_accessory)./reynolds_number + Kinf(vessel_accessory).*(1+(1./dapostrophie));
end

function dapostrophie = dapostrophie(vessel_diameter)
    % dapostrophie = 39.37 .* vessel_diameter;
    dapostrophie = vessel_diameter ./ 25.4;
end

function Kexp = Kexp(vessel_diameter_in, vessel_diameter_out)
    Kexp = (1-(vessel_diameter_in ./ vessel_diameter_out).^2).^2;
end

function Kcon = Kcon(vessel_diameter_in, vessel_diameter_out)
    Kcon = 0.5.*(1-(vessel_diameter_out ./ vessel_diameter_in).^2).^2;
end

function area = cross_sectional_area(vessel_diameter)
    area = 0.25 .* (pi .* (vessel_diameter.^2));
end

function reynolds_number = reynolds_number(vessel_viscosity, velocity_9_in, vessel_diameter_in, blood_density)
    reynolds_number = (blood_density.*velocity_9_in.*vessel_diameter_in) ./ vessel_viscosity;
end

function vessel_viscosity = vessel_viscosity(apparent_vessel_visc, Sj, H)
    vessel_viscosity = 1 + ((apparent_vessel_visc - 1) .* (((1-H).^Sj)./(0.55.^Sj)) );
end

function apparent_vessel_viscosity = apparent_vessel_viscosity(vessel_diameter)
    term1 = 2.20 .* exp(-1.3 .* vessel_diameter .* 10^-6);
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
    vol_flow = partion_coeff .* flow_coming_in;
end

function PC = partion_coeff_9(cross_sec_9_out, cross_sec_10_out)
    PC = cross_sec_9_out ./ (cross_sec_9_out + cross_sec_10_out);
end
