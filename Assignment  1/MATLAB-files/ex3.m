%% Setup
clear all, close all
%% constants
    m = 100;
    S = 14;
    C_d0 = 0.034;
    K = 0.07;
    g = 9.81;
    rho = 1.13;
    C_l = 0.7; %0.696932;
    c = [m S C_d0 K g rho C_l];

%% Initial values

tspan = linspace(0,100,2000);

stalling_speed = sqrt((g*2*m)/(C_l*S*rho))

n0 = 4;
x0 = linspace(0,0,n0);
h0 = linspace(10,10,n0);
v0 = linspace(1,stalling_speed*1.5,n0);
gamma0 = linspace(0,pi*1/4,n0);
t = zeros(length(tspan),n0);
y = zeros(length(tspan),n0,4);

opts = odeset('Events',@heightZeroEventsFcn);

%% ODE solver
for i = 1:n0
    y0 = [x0(i) h0(i) v0(i) gamma0(i)];
    [t_temp, y_temp, te_temp, ye_temp, ie_temp] = ode45(@(t,y) odefcn(y,c), tspan, y0, opts);
    %size(t_temp)
    t(1:length(t_temp),i) = t_temp;
    t(length(t_temp):end,i) = t_temp(end);
    y(1:length(y_temp),i,:) = y_temp;
    y(length(y_temp):end,i,:) = repmat(y_temp(end,:),length(tspan)- length(y_temp)+1,1);
    lables(i) = sprintf("$v_0 = %.1f m/s, \\gamma_0 = %.1f$", round(v0(i),1), round(rad2deg(gamma0(i)),1));
end

%% Plotting
figure;
hold on
plot(y(:,:,1),y(:,:,2))
xlabel('x');
ylabel('h');
ylim([0, inf])
legend(lables,'Interpreter','latex')
hold off

figure;
hold on
plot(t, y(:,:,3))
xlabel('$t$', 'Interpreter','latex');
ylabel('$\dot{v}$', 'Interpreter','latex');
legend(lables,'Interpreter','latex')

figure;
hold on
plot(t,y(:,:,3).*sin(y(:,:,4)))
xlabel('$t$', 'Interpreter','latex');
ylabel('$\dot{h}$', 'Interpreter','latex');
legend(lables,'Interpreter','latex')

figure;
hold on
plot(t,y(:,:,3).*cos(y(:,:,4)))
xlabel('$t$', 'Interpreter','latex');
ylabel('$\dot{x}$', 'Interpreter','latex');
legend(lables,'Interpreter','latex')

%% ODE function
function ret = odefcn(y, c)
    m = c(1);
    S = c(2);
    C_d0 = c(3);
    K = c(4);
    g = c(5);
    rho = c(6);
    C_l = c(7);
    x_dot = y(3)*cos(y(4));
    h_dot = y(3)*sin(y(4));
    v_dot = -(S*rho)/(2*m)*(C_d0 + K*C_l^2)*y(3)^2-g*sin(y(4));
    gamma_dot = 1/(2*m)*C_l*S*rho*y(3)-g/y(3)*cos(y(4));
    ret = [x_dot; h_dot; v_dot; gamma_dot];
end
%Stopping function
function [position,isterminal,direction] = heightZeroEventsFcn(t,y)
  position = y(2); % The value that we want to be zero
  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
end