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
u_0 = 35.7; % fuel injection [kg / s] <= Value changed in step block 
z1_0 = 1.0; % setting of valve z1 <= Value changed in step block 
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
fig_strings = ["step" "sin"] +"_extreme";

%% Battery

arr = [1 10 50 100 1000];
P_arr = P*0.12;
I_arr = I*0.0225;
D_arr = D *0.6;
P_2_arr = P_2*0.1;%*0.5;%arr
I_2_arr = I_2*0.3;%*arr
D_2_arr = D_2;%*arr
k_in_arr = k_in*20;%arr;	% amplification of battery charge control
k_out_arr = k_out*20;%arr;	% amplification of battery discharge control


step_response = 12.5;

sin_slow_f = 0.005;
sin_slow_a = 0.30*fkul_0;

sin_fast_f = 0.1;
sin_fast_a = sin_slow_a*(sin_slow_f/sin_fast_f);

use_step_arr = [1 0];

for use_step = use_step_arr
    for P = P_arr
        for I = I_arr
            for D = D_arr
                for P_2 = P_2_arr
                    for I_2 = I_2_arr
                        for D_2 = D_2_arr
                            for k_in = k_in_arr
                                for k_out = k_out_arr

%try
    fig_name = fig_strings(find(use_step == use_step_arr));
    %% Running the simulations inside of the Matlab-script
    Simulation_Time = 10000;
    SimOut = sim("voima_10.slx",Simulation_Time);
    
    %% figures;
    % u fig
    %figure('units','normalized','outerposition',[0 0 1 1]);
    figure();
%     sgtitle(sprintf("$P = %g, I = %g, D = %g \nP_2 = %g, I_2 = %g, D_2 = %g \nk_{in} =%g, k_{out} = %g $", ...
%                     arr(P == P_arr), ...
%                     arr(I == I_arr), ...
%                     arr(D == D_arr), ...
%                     arr(P_2 == P_2_arr), ...
%                     arr(I_2 == I_2_arr), ...
%                     arr(D_2 == D_2_arr), ...
%                     arr(k_in == k_in_arr), ...
%                     arr(k_out == k_out_arr)), ...
%         "Interpreter","latex","FontSize",13)
    %subplot(2,2,1)
    hold on
    plot(SimOut.time, SimOut.u,"LineWidth",1.5, "DisplayName", "$u$")
    ax = gca;
    ax.FontSize = 11;
    grid on
    %ylim([32 39])
    xlabel("$t$ [s]", "Interpreter","latex","FontSize",13);
    ylabel("$u$ [kg/s]", "Interpreter","latex","FontSize",13);
    legend("Interpreter","latex","FontSize",13)
    saveas(gcf, sprintf("Plots\\u_full_model_%s.png", fig_name))
    hold off

    % pkp fig  
    figure()
    %subplot(2,2,2)
    hold on
    plot(SimOut.time, SimOut.pkp,"LineWidth",1.5, "DisplayName", "$p_{kp}$")
    yline(1.02*pkp_0, "--r", "LineWidth",1.5, "DisplayName","$\pm2\%$")
    yline(0.98*pkp_0, "--r", "LineWidth",1.5, 'HandleVisibility','off')
    yline(pkp_0, "--b", "LineWidth",1.5, "DisplayName","$p_{kp0}$")
    ax = gca;
    ax.FontSize = 11;
    grid on
    ylim([0.97*pkp_0 1.03*pkp_0])
    xlabel("$t$ [s]", "Interpreter","latex","FontSize",13);
    ylabel("$p_{kp}$ [bar]", "Interpreter","latex","FontSize",13);
    legend("Interpreter","latex","FontSize",13)
    saveas(gcf, sprintf("Plots\\pkp_full_model_%s.png", fig_name))

    % pa fig  
    figure()
    %subplot(2,2,3)
    hold on
    plot(SimOut.time, SimOut.pa,"LineWidth",1.5, "DisplayName", "$p_{a}$")
    ax = gca;
    ax.FontSize = 11;
    grid on
    %ylim([0.97*pkp_0 1.03*pkp_0])
    xlabel("$t$ [s]", "Interpreter","latex","FontSize",13);
    ylabel("$p_{a}$ [bar]", "Interpreter","latex","FontSize",13);
    legend("Interpreter","latex","FontSize",13)
    saveas(gcf, sprintf("Plots\\pa_full_model_%s.png", fig_name))

    % consumption fig  
    figure()
    %subplot(2,2,4)
    hold on
    plot(SimOut.time, SimOut.c,"LineWidth",1.5, "DisplayName", "consumption")
    ax = gca;
    ax.FontSize = 11;
    grid on
    %ylim([0.97*pkp_0 1.03*pkp_0])
    xlabel("$t$ [s]", "Interpreter","latex","FontSize",13);
    ylabel("consumption", "Interpreter","latex","FontSize",13);
    legend("Interpreter","latex","FontSize",13)
    saveas(gcf, sprintf("Plots\\consumption_full_model_%s.png", fig_name))

                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
