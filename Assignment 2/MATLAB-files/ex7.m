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

    %% Running the simulations inside of the Matlab-script
    Simulation_Time = 1000;
    SimOut = sim("voima_7.slx",Simulation_Time);
    
    %% figures;
    % u fig
%     figure();
%     plot(SimOut.time, SimOut.u,"LineWidth",1.5)
%     ax = gca;
%     ax.FontSize = 11;
%     grid on
%     ylim([35.175 37])
%     xlabel("$t$ [s]", "Interpreter","latex","FontSize",13);
%     ylabel("$u$ [kg/s]", "Interpreter","latex","FontSize",13);
%     labels = ["$u$"];
%     legend(labels,"Interpreter","latex","FontSize",13)
%     saveas(gcf, sprintf("Plots\\u_tuning_%s.png", fig_name))
%     
%     % pkp fig
%     figure();
%     hold on
%     plot(SimOut.time, SimOut.pkp,"LineWidth",1.5)
%     yline([0.98*pkp_0 1.02*pkp_0], "--r", "LineWidth",1.5)
%     %
%     if(i==1)
%         dy=diff(SimOut.pkp)./diff(SimOut.time);
%         k=find(dy == max(dy));
%         tang=(SimOut.time-SimOut.time(k))*dy(k)+SimOut.pkp(k);
%         plot(SimOut.time,tang)
%         scatter(SimOut.time(k),SimOut.pkp(k))
%     end
%     hold off
%     ax = gca;
%     ax.FontSize = 11;
%     grid on
%     ylim([-inf, 93])
%     xlabel("$t$ [s]", "Interpreter","latex","FontSize",13);
%     ylabel("$p_{kp}$ [bar]", "Interpreter","latex","FontSize",13);
%     labels = ["$p_{kp}$" "$\pm2\%$"]; %Måste fixas så att båda linjerna är under samma 
%     if(i == 1)
%         labels = [labels "Tanget at steepst point" "Steepest point"];
%     end
%     legend(labels,"Interpreter","latex","FontSize",13)
%     saveas(gcf, sprintf("Plots\\pkp_tuning_%s.png", fig_name))
%     
    figure()
    hold on
    plot(SimOut.time, SimOut.pvp,"LineWidth",1.5, "DisplayName", "$p_{vp}$")
    yline(1.1*pvp_0, "--r", "LineWidth",1.5, "DisplayName","$\pm2\%$")
    yline(0.9*pvp_0, "--r", "LineWidth",1.5, 'HandleVisibility','off')
    ax = gca;
    ax.FontSize = 11;
    grid on
    xlabel("$t$ [s]", "Interpreter","latex","FontSize",13);
    ylabel("$p$ [bar]", "Interpreter","latex","FontSize",13);
    legend("Interpreter","latex","FontSize",13)
    %saveas(gcf, "Plots\\pvp_stepresponse.png")


