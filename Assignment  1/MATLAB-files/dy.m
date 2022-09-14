function Xdot = dy(X, cl)

% MS-E2133 System analysis laboratory II, assignment 1
% State equations of the glider

% parameters
m = 100;	    % mass
cd0 = 0.034;    % zero drag coefficient
k = 0.07;	    % induced drag coefficient
s = 14;		    % reference area
g = 9.809;	    % gravitational acceleration
rho = 1.13;	    % air density

% state variables
x = X(1);       % x-coordinate
h = X(2);       % altitude
v = X(3);       % velocity
gamma = X(4);   % flight path angle

% state equations
xdot = ??
hdot = ??
vdot = ??
gammadot = ??

Xdot = [xdot; hdot; vdot; gammadot];
