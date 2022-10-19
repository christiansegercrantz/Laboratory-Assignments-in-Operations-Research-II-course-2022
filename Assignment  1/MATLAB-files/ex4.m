%% Setup
clear all, close all
%% constants
m = 100;
S = 14;
C_D0 = 0.034;
K = 0.07;
g = 9.81;
rho = 1.13;
%C_l = 1.4;
%c = [m S C_d0 K g rho C_l];

%% Values
C_L = linspace(-1.4,1.4,100);
C_D = C_D0 + K*C_L.^2;
tangent = 2*K*sqrt(C_D0/K)*C_L;
x0 = sqrt(C_D0/K);

plot(C_L,C_D,'LineWidth',1.5)
ylim([0,0.2])
grid on
hold on
plot(C_L,tangent,'LineWidth',1.5)
scatter(x0, 2*K*sqrt(C_D0/K)*x0,'LineWidth',1.5)
ax = gca;
ax.FontSize = 11;
lables = ["$C_{D_0}+KC_L^2$", "$C_L$", strcat("$y=",num2str(2*K*sqrt(C_D0/K)*x0),", x=", num2str(x0), "$")]
xlabel('$C_L$', 'Interpreter','latex','FontSize',13);
ylabel('$C_D$', 'Interpreter','latex','FontSize',13);
legend(lables,'Interpreter','latex','FontSize',13)
saveas(gcf, "Plots\optimal_angle.png")