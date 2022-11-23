% ==================================================
% MS-E2133 - Systems analysis laboratory II
% Matlab-script: Simulation and control of a thermal power plant
% ==================================================
% 
% Run this script to initialize the values needed in model.mdl
% 
% See example in lines 50-54 for running the model inside of Matlab

%% Parameters and initial values

T1 = 30.0; % time constant of the steam production in the boiler [s]
Ts = 15.0; % storage capacity of the boiler [kg/bar]

% Conversion factors:
k0 = 1.0; % conversion factor from the fuel injection to the steam production
k1 = 2.0;	% flow resistance from the boiler to the steam stock [kg/(s*bar)]
k2 = 1/30;	% storage capacity of the steam stock [kg/bar]
k3 = 0.4;	% flow resistance to the turbine [kg/(s*bar)]
k4 = 1/300;	% storage capacity of the counter pressure stock [kg/bar]
k5 = 0.6;	% flow resistance into the steam battery [kg/(s*bar)]
k6 = 0.6;	% flow resistance out of the steam battery [kg/(s*bar)]

% Steam battery:
h1 = 8.3; % parameter of the steam battery
h2 = 21.0; % parameter of the steam battery [bar]
pcp = 13.0;	% battery charging pressure
k_in = 25.5;	% amplification of battery charge control
k_out = 210.0;	% amplification of battery discharge control


% Equilibrium values:
% u_0 = 35.7; % fuel injection [kg / s] <= Value changed in step block 
% z1_0 = 1.0; % setting of valve z1 <= Value changed in step block 
fp_0   = 35.7;	% steam generation in the boiler [kg/s]
pk_0   = 91.0;	% pressure of the boiler [bar]
pkp_0  = 89.25;	% pressure of the high pressure stock [bar]


fin_0  = 0.0;	% flow into the battery [kg/s]
fout_0 = 0.0;	% flow out of the battery [kg/s]
pa_0   = 10.0;	% pressure of the battery [bar]
ma_0   = 150000.0; % mass of the water in the battery [kg]

fkul_0 = 35.7; % flow to consumption [kg/s]
pvp_0  = 3.0;	% pressure of the counter pressure stock [bar]


%% Running the simulations inside of the Matlab-script
Simulation_Time = 1000;
SimOut = sim('voima2',Simulation_Time);

%% figures;
% u fig
figure();
plot(SimOut.time, SimOut.u,'LineWidth',1.5)
ax = gca;
ax.FontSize = 11;
grid on
%title('control')
xlabel('$t$ [s]', 'Interpreter','latex','FontSize',13);
ylabel('$u$ [kg/s]', 'Interpreter','latex','FontSize',13);

% pkp fig
figure();
plot(SimOut.time, SimOut.pkp,'LineWidth',1.5)
ax = gca;
ax.FontSize = 11;
grid on
xlabel('$t$ [s]', 'Interpreter','latex','FontSize',13);
ylabel('$p_{kp} [bar]$', 'Interpreter','latex','FontSize',13);
saveas(gcf, "Plots\\pkp_openvalve.png")

% fp fig
figure();
plot(SimOut.time, SimOut.fp,'LineWidth',1.5)
ax = gca;
ax.FontSize = 11;
grid on
xlabel('$t$ [s]', 'Interpreter','latex','FontSize',13);
ylabel('$f_{p} [kg/s]$', 'Interpreter','latex','FontSize',13);
saveas(gcf, "Plots\\fp_openvalve.png")

% fg fig
figure();
plot(SimOut.time, SimOut.fg,'LineWidth',1.5)
ax = gca;
ax.FontSize = 11;
grid on
xlabel('$t$ [s]', 'Interpreter','latex','FontSize',13);
ylabel('$f_{g} [kg/s]$', 'Interpreter','latex','FontSize',13);
saveas(gcf, "Plots\\fg_openvalve.png")

% z1 fig
figure();
plot(SimOut.time, SimOut.z1,'LineWidth',1.5)
ax = gca;
ax.FontSize = 11;
grid on
xlabel('$t$ [s]', 'Interpreter','latex','FontSize',13);
ylabel('$z_{1}$', 'Interpreter','latex','FontSize',13);
saveas(gcf, "Plots\\z1_openvalve.png")
