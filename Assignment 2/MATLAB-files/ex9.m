% ==================================================
% MS-E2133 - Systems analysis laboratory II
% Matlab-script: Simulation and control of a thermal power plant
% ==================================================
% 
% Run this script to initialize the values needed in model.mdl
% 
% See example in lines 50-54 for running the model inside of Matlab
%%
close all  


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


%% PID values
a = pkp_0-88.9958;
tau = 16;
D =  1.2/a*0.5*tau;
I = (1.2/a)/(2*tau);
P = 1.2/a;

%% PID2 values
P_2 = 100;
I_2 = 0.1; 
D_2 = 0;
fig_strings = ["equal" "k_in_small" "k_out_small" "k_in_large" "k_out_large" "k_out_huge"];

k_in_arr  =  k_in*[1 0.1 1   10 1 1  ];	% amplification of battery charge control
k_out_arr = k_out*[1 1   0.1 1 10 250];	% amplification of battery discharge control

for i = 1:length(k_in_arr)
    k_in = k_in_arr(i)
    k_out = k_out_arr(i)
    fig_name = fig_strings(i);
    %% Running the simulations inside of the Matlab-script
    if i <= 5
        Simulation_Time = 1000;
    else
        Simulation_Time = 35000;
    end
    
    SimOut = sim("voima_9.slx",Simulation_Time);
    
    %% figures;
    % u fig
    figure()
    %figure('units','normalized','outerposition',[0 0 1 1]);
    %subplot(2,1,1)
    hold on
    plot(SimOut.time, SimOut.u,"LineWidth",1.5, "DisplayName", sprintf("$u, k_{in} = %g, k_{out} = %g$", k_in, k_out))
    %hold off
    ax = gca;
    ax.FontSize = 11;
    grid on
    ylim([35.175 41])
    xlabel("$t$ [s]", "Interpreter","latex","FontSize",13);
    ylabel("$u$ [kg/s]", "Interpreter","latex","FontSize",13);
    legend("Interpreter","latex","FontSize",13)
    saveas(gcf, sprintf("Plots\\u_steam_battery_%s.png", fig_name))

    % pa fig  
    figure()
    %subplot(2,1,2)
    hold on
    plot(SimOut.time, SimOut.pa,"LineWidth",1.5, "DisplayName", sprintf("$p_{a}, k_{in} = %g, k_{out} = %g$", k_in, k_out))
    ax = gca;
    ax.FontSize = 11;
    grid on
    if i <= 4
        ylim([9.95 10.1]);
    end
    xlabel("$t$ [s]", "Interpreter","latex","FontSize",13);
    ylabel("$p_a$ [bar]", "Interpreter","latex","FontSize",13);
    legend("Interpreter","latex","FontSize",13)
    saveas(gcf, sprintf("Plots\\pa_steam_battery_%s.png", fig_name))
end
