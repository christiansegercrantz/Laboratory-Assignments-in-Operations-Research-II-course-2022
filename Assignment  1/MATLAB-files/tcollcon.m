function [g, geq] = tcollcon(X, dp, sc)

% MS-E2133 System analysis laboratory II, assignment 1
% Optimal flight with a glider
% A function required by the Matlab-script 'fmincon', that returns the
% values of the constraints at a given point.

g = [];    % inequality constraints

tf = X(end) * sc(5);    % terminal time (s)
dt = tf / (dp - 1);	% length of the discretization step (s)

% CONSTRAINTS:
%=============

%1 form the state matrix
Xm = zeros(4, dp);
for i = 1:4 		    % state variables (x, h, vx, vy)
    Xm(i, :) = X((i-1)*dp+1:i*dp) * sc(i);
end

%2 form the controls
n = X(4*dp+1:end);

%3 calculate xdot st the discretization points
Xdotm = zeros(4, dp);
for j = 1:dp
  Xdotm(:, j) = tdy(Xm(:,j), n(j)); % state and control
end

%4 calculate the piecewise polynomial and its derivative at the centres of
%the intervals
Xmidm = zeros(4, dp-1);
Xdotmidm = zeros(4, dp-1);
for i = 1:4
    Xmidm(i, :) = (Xm(i, 1:dp-1) + Xm(i, 2:dp))/2 +...
        dt*(Xdotm(i, (1:dp-1)) - Xdotm(i, 2:dp))/8;
    Xdotmidm(i, :) = -3*(Xm(i, (1:dp-1)) - Xm(i, (2:dp)))/(2*dt) -...
        (Xdotm(i, (1:dp-1)) + Xdotm(i, (2:dp)))/4;
end         

%5 calculate the state equations of the system at the center points
systmidm = zeros(4, dp - 1);
for j = 1:dp-1
  systmidm(:, j) = tdy(Xmidm(:, j), (n(j+1) + n(j))/2);
end

%6 finally obtain the collocation constraints
%  (control constraints are handled as simple bounds)
csc=[10 10 1 1];  %scaling
geq = zeros(1, 4*(dp-1)); % equality constraints (systmidm=Xdotmidm)
for i = 1:4									
    geq((i-1)*(dp-1)+1:i*(dp-1)) =...
        (systmidm(i, 1:dp-1) - Xdotmidm(i, 1:dp-1)) ./ csc(i);
end
