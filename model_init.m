clc; clear all; close all;

%% Vessel constants
vessel.L   = 4.8;
vessel.B   = 0.78;
vessel.T   = 0.24;
vessel.Cb  = 0.3;
vessel.R66 = 0.25*vessel.L;
vessel.cg  = [0 0 0]';
vessel.m   = 280;
vessel.T_surge = 1;
vessel.T_sway  = 2;
vessel.T_yaw   = 1;

% Rudder Model
vessel.H = 0.42;                  % Rudder plate height.
vessel.W = 0.25;                  % Rudder plate width.
vessel.rho_water = 1025;          % Water density
vessel.A_R = vessel.H*vessel.W;   % Rudder area, [m^2]

vessel.rudder_saturation = 45*(pi/180);

% Fossen rudder forces equations
vessel.lambda    = vessel.H^2 / vessel.A_R;
vessel.C_N       = 6.13 / (vessel.lambda + 2.25);
vessel.t_R       = 1 - (0.28*vessel.Cb+0.55);
vessel.a_H       = 0.2;                      % From graphs in Fossen (2021)
vessel.x_prime_H = -1.8;               % From graphs in Fossen (2021)
vessel.x_H       = vessel.x_prime_H / vessel.L;
vessel.x_R       = -0.5 * vessel.L;

% Wind coefficients
vessel.rho_air = 1.2;    % Density, air
vessel.L_oa    = 5;         % Length overall
vessel.A_fw    = 0.7 * 0.24; % Frontal projected area
vessel.A_lw    = vessel.L_oa * 0.24; % Lateral projected area
vessel.c_x     = 0.5;       % [0.50,0.90]
vessel.c_y     = 0.7;       % [0.70,0.95]
vessel.c_n     = 0.05;      % [0.05,0.20]



%% Compute M and D
Iz = 1/12*vessel.m*(vessel.L^2 + vessel.B^2);

MRB = [      vessel.m                 0              -vessel.m*vessel.cg(2);           % rigid-body inertia matrix
                0                  vessel.m           vessel.m*vessel.cg(1);                
       -vessel.m*vessel.cg(2)  vessel.m*vessel.cg(1)           Iz];


vessel.MRB = MRB;
% Nominal added mass using Clark83
% Nondimenisonal hydrodynamic derivatives in surge
U = 0.5;
Xudot = -0.1 * vessel.m / (0.5 * vessel.rho_water * vessel.L^3); 
Xu = -((vessel.m-Xudot) / vessel.T_surge) / (0.5 * vessel.rho_water  * vessel.L^2 * U);  
 
% Nondimenisonal hydrodynamic derivatives in sway and yaw from Clarke et al. (1983)
S = pi * (vessel.T/vessel.L)^2;                 % scale factor 

Yvdot = -S * ( 1 + 0.16 * vessel.Cb * vessel.B/vessel.T - 5.1 * (vessel.B/vessel.L)^2 );
Yrdot = -S * ( 0.67 * vessel.B/vessel.L - 0.0033 * (vessel.B/vessel.T)^2 );
Nvdot = -S * ( 1.1 * vessel.B/vessel.L - 0.041 * (vessel.B/vessel.T) );
Nrdot = -S * ( 1/12 + 0.017 * vessel.Cb * (vessel.B/vessel.T) - 0.33 * (vessel.B/vessel.L) );
Yv = -S * ( 1 + 0.4 * vessel.Cb * (vessel.B/vessel.T) );
Yr = -S * ( -1/2 + 2.2 * (vessel.B/vessel.L) - 0.08 * (vessel.B/vessel.T) );
Nv = -S * ( 1/2 + 2.4 * (vessel.T/vessel.L) );
Nr = -S * ( 1/4 + 0.039 * (vessel.B/vessel.T) - 0.56 * (vessel.B/vessel.L) );
 
% Nondimenisonal hydrodynamic matrices 
MA_prime = [ -Xudot   0        0
               0      -Yvdot   -Yrdot
               0      -Nvdot   -Nrdot ];

 
% Dimensional model (Fossen 2021, Appendix D)   
T    = diag([1 1 1/vessel.L]);
Tinv = diag([1 1 vessel.L]);

MA = (0.5 * vessel.rho_water * vessel.L^3) * Tinv^2 * (vessel.T * MA_prime * Tinv);
vessel.MA  = (MA + MA')/2;
vessel.M = vessel.MRB + vessel.MA;
vessel.D = diag(diag(vessel.M) ./ [vessel.T_surge; vessel.T_sway; vessel.T_yaw]);

%% Speed model
speed_model = [ 0.116392998053662;...
                0.214487083945715;...
                0.0880678632611925;...
               -0.00635496887217675;...
                0.0937464223577265;...
                0.238364678400396];


%% Environmental parameters
% Waves
Hs = 3;
Tp = 10;
wave_dir = 90*(pi/180);

% Wind
wind_speed     = 5;
wind_direction = 90*(pi/180); % Wind direction going to

% Current
current_speed = 0.0;
current_direction = 90*(pi/180);


%% Simulation parameters
sim_time = 50;
nu0  = [1;0;0];
eta0 = zeros(3,1);

logging_frequency = 10; % Hz

sim_res = sim('model_3DOF');


%% Plots

t = sim_res.desired_heading.time;
desired_heading = sim_res.desired_heading.signals.values;

N = squeeze(sim_res.eta.signals.values(1,1,:));
E = squeeze(sim_res.eta.signals.values(2,1,:));
psi = squeeze(sim_res.eta.signals.values(3,1,:));
u = squeeze(sim_res.nu.signals.values(1,1,:));
v = squeeze(sim_res.nu.signals.values(2,1,:));
r = squeeze(sim_res.nu.signals.values(3,1,:));

sog = sim_res.sog.signals.values;
cog = sim_res.cog.signals.values;


figure
plot(t,desired_heading); grid on; hold on
plot(t,psi)
plot(t,cog)
xlabel('Time [s]')
ylabel('Heading [rad]')
legend('Desired heading','True heading','COG')

figure
subplot(311)
plot(t,u); grid on
xlabel('Time [s]')
ylabel('Surge velocity [m/s]')
subplot(312)
plot(t,v); grid on
xlabel('Time [s]')
ylabel('Sway velocity [m/s]')
subplot(313)
plot(t,r); grid on
xlabel('Time [s]')
ylabel('Yaw angular velocity [rad/s]')

figure
plot(E,N); grid on
xlabel('East [m]')
ylabel('North [m]')
