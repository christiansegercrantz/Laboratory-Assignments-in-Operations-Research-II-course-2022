function Xdot = tdy(X, cl)

% MS-E2133 System analysis laboratory II, assignment 1
% State equations of the glider

% parameters
m = 100;	    % mass
cd0 = 0.034;	% zero drag coefficient
k = 0.07;	    % induced drag coefficient
s = 14;	    	% reference area
g = 9.809;	    % gravitational acceleration
rho = 1.13;	    % air density
umax = 2.5;	    % maximum air flow upwards in thermal
R = 100;		% constant (m)
x0 = 150;		% constant (m)

% state variables
x = X(1);       % x-coordinate
h = X(2);       % altitude
vx = X(3);      % x-velocity
vh = X(4);      % h-velocity

ua = umax*exp(-((x - x0)/R)^2)*(1 - ((x - x0)/R)^2); % air flow upward in thermal (m/s)
vr = sqrt(vx^2 + (vh - ua)^2);                       % velocity relative to air
L = 0.5*cl*rho*s*vr^2;			                     % lift
cd = cd0 + k*cl^2;                                   % polar
D = 0.5*cd*rho*s*vr^2;			                     % resistance force
singa = (vh - ua)/vr;
cosga = vx/vr;

% state equations
xdot = vx;
hdot = vh;
vxdot = (-L*singa-D*cosga)/m;
vhdot = (L*cosga-D*singa-m*g)/m;

Xdot = [xdot; hdot; vxdot; vhdot];
