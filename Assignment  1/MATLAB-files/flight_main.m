
% ==================================================
% MS-E2133 System analysis laboratory II
% Optimal flight with a glider
% ==================================================

% State variables are x, h, v and gamma

clear
close all

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
    figure(1)
    clf
    
    subplot(221)
    plot(a(1:iter)*sc(1) ,a(iter+1:2*iter)*sc(2), 'b-+')
    title('altitude vs. position')
    grid on
    xlabel('{\itx}, m')
    ylabel('{\ith}, m')
    
    subplot(222)
    plot(linspace(0, a(end), iter)*sc(5), a(2*iter+1:3*iter)*sc(3), 'b-+')
    grid on
    title('velocity vs. time')
    xlabel('{\itt}, s')
    ylabel('{\itv}, m/s')
    
    subplot(223)
    ka = (a(4*iter+1:5*iter-1) + a(4*iter+2:5*iter)) / 2;
    plot(a(end)*linspace(1/(iter-1)/2, 1-1/(iter-1)/2, iter-1)*sc(5), ka, 'b-+');
    grid on
    title('control')
    xlabel('{\itt}, s')
    ylabel('{\itC_L}')
    
    subplot(224)
    plot(linspace(0, a(end), iter)*sc(5), a(3*iter+1:4*iter), 'b-+')
    grid on
    title('flight path angle vs. time')
    xlabel('{\itt}, s')
    ylabel('{\it\gamma}, rad')
    
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
