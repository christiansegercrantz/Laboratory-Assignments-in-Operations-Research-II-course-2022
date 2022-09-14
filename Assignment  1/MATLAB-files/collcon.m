function [g, geq] = collcon(X, dp, sc)

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
for ii = 1:4 		    % state variables (x, h, v, gamma)
    Xm(ii, :) = X((ii-1)*dp+1:ii*dp) * sc(ii);
end

%2 form the controls
n = X(4*dp+1:end);

%3 calculate xdot st the discretization points
Xdotm = zeros(4, dp);
for j = 1:dp
  Xdotm(:, j) = dy(Xm(:,j), n(j)); % state and control
end

%4 calculate the piecewise polynomial and its derivative at the centres of
%the intervals
Xmidm = zeros(4, dp - 1);
Xdotmidm = zeros(4, dp - 1);
for ii = 1:4
    Xmidm(ii, :) = (Xm(ii, 1:dp-1) + Xm(ii, 2:dp))/2 +...
        dt*(Xdotm(ii, (1:dp-1)) - Xdotm(ii, 2:dp))/8;
    Xdotmidm(ii, :) = -3*(Xm(ii, (1:dp-1)) - Xm(ii, (2:dp)))/(2*dt) -...
        (Xdotm(ii, (1:dp-1)) + Xdotm(ii, (2:dp)))/4;
end         

%5 calculate the state equations of the system at the center points
systmidm = zeros(4, dp - 1);
for j = 1:dp-1
  systmidm(:, j) = dy(Xmidm(:, j), (n(j+1) + n(j))/2);
end

%6 finally obtain the collocation constraints
%  (control constraints are handled as simple bounds)
csc = [1, 1, 1, 1];  	%Possibility for scaling.
					%Scaling is not applied here, so all values are 1.

geq = zeros(1, 4*(dp-1)); % equality constraints (systmidm=Xdotmidm)
for ii = 1:4									
    geq((ii-1)*(dp-1)+1:ii*(dp-1)) =...
        (systmidm(ii, 1:dp-1) - Xdotmidm(ii, 1:dp-1)) ./ csc(ii);
end
