clear
close all

A = dlmread('veratk_2.txt');
t_compare = A(1,:);
x_compare = A(2,:);
h_compare = A(3,:);
v_compare = A(4,:);
gamma_compare = A(5,:);
px_compare = A(6,:);
ph_compare = A(7,:);
pv_compare = A(8,:);
pgamma_compare = A(9,:);

dp = 5;     % number of discretization points at the start
step = 3;   % increment of discretization points
dpend = 20; % number of discretization points at the end

sc = [150, 50, 20, 1, 30];    %scaling (x, h, v, gamma, t)
%sc = [1, 1, 1, 1, 1];

% initial and terminal conditions
x_0 = 0;
h_0 = 50;
h_f = 40;
v_0 = 13;
gamma_0 = 0;

% settings of the optimization algorithm
options = optimset(...
    'LargeScale', 'off',...
    'Display','iter',...
    'TolFun', 1d-3,...
    'TolX', 1d-3,...
    'TolCon', 1d-4,...
    'MaxFunEvals', 10000,...
    'Algorithm', 'sqp');

% initial guesses
x = linspace(x_0, 300, dp);
h = linspace(h_0, h_f, dp);
v = v_0 * ones(1, dp);
gamma = gamma_0 * ones(1, dp);
cl = ones(1,dp);
X = [x/sc(1), h/sc(2), v/sc(3), gamma/sc(4), cl, 20/sc(5)];


% continuation loop
% =================

for iter = dp:step:dpend
    
    % upper bounds for variables
    ub_x = 1000 / sc(1) * ones(1, iter);
    ub_x(1) = x_0 / sc(1);
    
    ub_h = 90 / sc(2) * ones(1, iter);
    ub_h(1) = h_0 / sc(2);
    ub_h(iter) = h_f / sc(2);
    
    ub_v = 40 / sc(3) * ones(1, iter);
    ub_v(1) = v_0 / sc(3);
    
    ub_gamma = 1.5 / sc(4) * ones(1, iter);
    ub_gamma(1) = gamma_0 / sc(4);
    
    ub_cl = 1.4 * ones(1, iter);
    
    ub_t = 50 / sc(5);
    
    ub = [ub_x, ub_h, ub_v, ub_gamma, ub_cl, ub_t];
    
    % lower bounds for variables
    lb_x = zeros(1, iter);
    lb_x(1) = ub_x(1);
    
    lb_h = zeros(1, iter);
    lb_h(1) = ub_h(1);
    lb_h(iter) = ub_h(iter);
    
    lb_v = 5 / sc(3) * ones(1, iter);
    lb_v(1) = ub_v(1);
    lb_v(iter) = 10 / sc(3);
    
    lb_gmmaa = -1.5 / sc(4) * ones(1, iter);
    lb_gmmaa(1) = ub_gamma(1);
    
    lb_cl = -1.4 * ones(1, iter);
    
    lb_t = 1 / sc(5);
    
    lb = [lb_x lb_h lb_v lb_gmmaa lb_cl lb_t];
    
    % optimization
    [a, xf, exitflag, output, lambda] = fmincon('objfun', X,...
        [], [], [], [], lb, ub, 'collcon', options, iter, sc);
    
    % interpolation for new initial guesses
    origtime = linspace(0, a(end), iter);
    itime = linspace(0, a(end), iter + step);
    X_x = interp1(origtime, a(1:iter), itime, 'linear');
    X_h = interp1(origtime, a(iter+1:2*iter), itime, 'linear');
    X_v = interp1(origtime, a(2*iter+1:3*iter), itime, 'spline');
    X_gamma = interp1(origtime, a(3*iter+1:4*iter), itime, 'spline');
    cl = interp1(origtime, a(4*iter+1:5*iter), itime, 'spline');
    X = [X_x, X_h, X_v, X_gamma, cl, a(end)];
    
    % plot the results
    %figure(1)
    %clf
    
    figure(1)%subplot(221)
    plot(a(1:iter)*sc(1) ,a(iter+1:2*iter)*sc(2), 'b-+')
    ax = gca;
    ax.FontSize = 11;
    %title('altitude vs. position')
    grid on
    xlabel('$x$ [m]', 'Interpreter','latex','FontSize',13);
    ylabel('$h$ [m]', 'Interpreter','latex','FontSize',13);
    
    figure(2)%subplot(222)
    plot(linspace(0, a(end), iter)*sc(5), a(2*iter+1:3*iter)*sc(3), 'b-+')
    ax = gca;
    ax.FontSize = 11;
    grid on
    %title('velocity vs. time')
    xlabel('$t$ [s]', 'Interpreter','latex','FontSize',13);
    ylabel('$v$ [m/s]', 'Interpreter','latex','FontSize',13);
    
    figure(3)%subplot(223)
    ka = (a(4*iter+1:5*iter-1) + a(4*iter+2:5*iter)) / 2;
    plot(a(end)*linspace(1/(iter-1)/2, 1-1/(iter-1)/2, iter-1)*sc(5), ka, 'b-+');
    ax = gca;
    ax.FontSize = 11;
    grid on
    %title('control')
    xlabel('$t$ [s]', 'Interpreter','latex','FontSize',13);
    ylabel('$C_L$', 'Interpreter','latex','FontSize',13);
    
    figure(4)%subplot(224)
    plot(linspace(0, a(end), iter)*sc(5), rad2deg(a(3*iter+1:4*iter)), 'b-+')
    grid on
    ax = gca;
    ax.FontSize = 11;
    %title('flight path angle vs. time')
    xlabel('$t$ [s]', 'Interpreter','latex','FontSize',13);
    ylabel('$\gamma$ [deg]', 'Interpreter','latex','FontSize',13);
    
    drawnow;
    
    fprintf('\nThe number of discretization points was %.0f\n', iter);
    fprintf('x(tf) = %.2f m\n', a(iter)*sc(1));
    fprintf('Distance travelled: x(tf) - x(0) = %.2f m\n', (a(iter)-a(1))*sc(1));
    fprintf('v(tf) = %.2f m/s\n', a(3*iter)*sc(3));
    fprintf('Terminal time tf = %.2f s\n\n\n', a(end)*sc(5));
end

% print Lagrange multipliers
lx = lambda.eqnonlin(1:dpend-1);
lh = lambda.eqnonlin(dpend:2*dpend-2);
lv = lambda.eqnonlin(2*dpend-1:3*dpend-3);
lg = lambda.eqnonlin(3*dpend-2:4*dpend-4);
L = [lx'; lh'; lv'; lg'];
fprintf('Lagrange multipliers: L=\n');
fprintf('%.4f ',lx);
fprintf('\n');
fprintf('%.4f ',lh);
fprintf('\n');
fprintf('%.4f ',lv);
fprintf('\n');
fprintf('%.4f ',lg);
fprintf('\n');

t = linspace(0, a(end), iter-1)*sc(5);
delta_t = (dpend-1)/(a(end)*sc(5));
lables = ["Dynamically obtained result", "Results to compare to"];

figure(1)
hold on
plot(x_compare, h_compare, "r--")
hold off
legend(lables,'Interpreter','latex','FontSize',13, 'Location','southwest')
saveas(figure(1), "Plots\dynamic_h_compare.png")

figure(2)
hold on
plot(t_compare, v_compare, "r--")
hold off
legend(lables,'Interpreter','latex','FontSize',13, 'Location','southwest')
saveas(figure(2), "Plots\dynamic_v_compare.png")

C_L_compare = pgamma_compare./(2*0.07*v_compare.*pv_compare);

figure(3)%subplot(223)
hold on
plot(t_compare, C_L_compare, "r--")
hold off
legend(lables,'Interpreter','latex','FontSize',13, 'Location','southwest')
saveas(figure(3), "Plots\dynamic_control_compare.png")

figure(4)
hold on
plot(t_compare, rad2deg(gamma_compare), "r--")
hold off
legend(lables,'Interpreter','latex','FontSize',13, 'Location','southwest')
saveas(figure(4), "Plots\dynamic_gamma_compare.png")

figure;
hold on
grid on
plot(t, 150*lx*3/2*delta_t, "-x")
plot(t_compare, px_compare, "-o")
ax = gca;
ax.FontSize = 11;
lables = ["$L_x$", "$p_x$"];
legend(lables,'Interpreter','latex','FontSize',13, 'Location','west')
saveas(gcf, "Plots\costates_compare_px.png")

figure;
hold on
grid on
plot(t, 150*lh*3/2*delta_t, "-x")
plot(t_compare,ph_compare, "-o")
ax = gca;
ax.FontSize = 11;
lables = ["$L_h$", "$p_h$"];
legend(lables,'Interpreter','latex','FontSize',13, 'Location','west')
saveas(gcf, "Plots\costates_compare_ph.png")

figure;
hold on
grid on
plot(t, 150*lv*3/2*delta_t, "-x")
plot(t_compare, pv_compare, "-o")
ax = gca;
ax.FontSize = 11;
lables = ["$L_v$", "$p_v$"];
legend(lables,'Interpreter','latex','FontSize',13, 'Location','northwest')
saveas(gcf, "Plots\costates_compare_pv.png")

figure;
hold on
grid on
plot(t, 150*lg*3/2*delta_t, "-x")
plot(t_compare, pgamma_compare, "-o")
ax = gca;
ax.FontSize = 11;
lables = ["$L_g$", "$p_g$"];
legend(lables,'Interpreter','latex','FontSize',13, 'Location','northwest')
saveas(gcf, "Plots\costates_compare_pg.png")

