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


    %% State controller values
    A = [-1/T1 0 0 0
         1/Ts -(k1*pk_0)/(sqrt(pk_0^2-pkp_0^2)*Ts) (k1*pk_0)/(sqrt(pk_0^2-pkp_0^2)*Ts) 0
         0 (k1*k2*pk_0)/(sqrt(pk_0^2-pkp_0^2)) k2*(-(k1*pk_0)/(sqrt(pk_0^2-pkp_0^2))-k3*z1_0) 0
         0 0 -1 0];
    
    B = [k0/T1
        0
        0
        0];
    
    %C = [0 0 1 0];
    
    %D = 0;
    
    Qconst = 1/(pkp_0*1.02)^2;
    Rconst = 1/(u_0*1.02)^2;
    
    n=4;
    Q11 = 1e-6;
    Q22 = 1;
    Q33 = Qconst/100;
    Q44 = 1e-4;
    
    R1 = 0.06;


    Q = [ Q11 0    0      0
          0   Q22  0      0
          0   0       Q33 0
          0   0       0      Q44];
    
    R = R1;

    [K, S, P] = lqr(A,B,Q,R);

    g1= K(1); g2 = K(2); g3 = K(3); g4 = K(4);

    %% PID values
    P = 5.66;
    I = 0.1;
    D = 37.77;
    
    %Setting
    use_pid_arr = [1 0];
    fig_strings = ["PID" "State"] ;

for use_pid = use_pid_arr
    fig_name = fig_strings(find(use_pid == use_pid_arr));
    %% Running the simulations inside of the Matlab-script
    Simulation_Time = 1000;
    SimOut = sim("voima_6.slx",Simulation_Time);
    
    %% figures;
    % u fig
    figure(1);
    hold on
    plot(SimOut.time, SimOut.u,"LineWidth",1.5, "DisplayName", sprintf("%s", fig_name))
    %plot(SimOut.time, SimOut.u,"LineWidth",1.5, "DisplayName", sprintf("$R = %g$", R))
    %hold off
    ax = gca;
    ax.FontSize = 11;
    grid on
    %ylim([35.175 37])
    xlabel("$t$ [s]", "Interpreter","latex","FontSize",13);
    ylabel("$u$ [kg/s]", "Interpreter","latex","FontSize",13);
    legend("Interpreter","latex","FontSize",13)
    
    
    % pkp fig
    figure(2);
    hold on
    if use_pid == 1
        yline(1.02*pkp_0, "--r", "LineWidth",1.5, "DisplayName","$\pm2\%$")
        yline(0.98*pkp_0, "--r", "LineWidth",1.5, 'HandleVisibility','off')
    end
    plot(SimOut.time, SimOut.pkp,"LineWidth",1.5, 'DisplayName', sprintf("%s", fig_name))
    %plot(SimOut.time, SimOut.pkp,"LineWidth",1.5, 'DisplayName', sprintf("$R = %g$", R))
    %hold off
    ax = gca;
    ax.FontSize = 11;
    grid on
    ylim([-inf, 93])
    xlabel("$t$ [s]", "Interpreter","latex","FontSize",13);
    ylabel("$p_{kp}$ [bar]", "Interpreter","latex","FontSize",13);
    legend("Interpreter","latex","FontSize",13)
    %saveas(gcf, sprintf("Plots\\pkp_statefeedback_integral_%s.png", int2str(i)))
    
%     % fp fig
%     figure();
%     hold on
%     plot(SimOut.time, SimOut.fp,"LineWidth",1.5, "DisplayName", sprintf("$Q_{44} = %g$", Q44(i)))
%     ax = gca;
%     ax.FontSize = 11;
%     grid on
%     xlabel("$t$ [s]", "Interpreter","latex","FontSize",13);
%     ylabel("$f_{p}$ [kg/s]", "Interpreter","latex","FontSize",13);
%     legend("Interpreter","latex","FontSize",13)
%     %saveas(gcf, sprintf("Plots\\fp_statefeedback_integral_%s.png", int2str(i)))
%     
%     % pk fig
%     figure()
%     hold on
%     plot(SimOut.time, SimOut.pk,"LineWidth",1.5, "DisplayName", sprintf("$Q_{44} = %g$", Q44(i)))
%     ax = gca;
%     ax.FontSize = 11;
%     grid on
%     xlabel("$t$ [s]", "Interpreter","latex","FontSize",13);
%     ylabel("$k_{k}$ [bar]", "Interpreter","latex","FontSize",13);
%     legend("Interpreter","latex","FontSize",13)
%     %saveas(gcf, sprintf("Plots\\pk_statefeedback_integral_%s.png", int2str(i)))

    % 
    % % z1 fig
    % fig_z1 = figure();
    % plot(SimOut.time, SimOut.z1,"LineWidth",1.5)
    % ax = gca;
    % ax.FontSize = 11;
    % grid on
    % xlabel("$t$ [s]", "Interpreter","latex","FontSize",13);
    % ylabel("$z_{1}$", "Interpreter","latex","FontSize",13);
end
saveas(figure(1), "Plots\\u_comparison.png")
saveas(figure(2), "Plots\\pkp_comparison.png")