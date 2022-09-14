function f = tobjfun(X, dp, sc)

% MS-E2133 System analysis laboratory II, assignment 1
% Optimal flight with a glider
% A function required by the Matlab-script 'fmincon', that returns the
% value of the function to be minimized at a given point.

% The function value is the x-coordinate at the end (fmincon minimizes)
f = -X(dp);