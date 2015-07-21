function [T, Y] = two_body_caller(pos, vel, err, bigmass, time)
% Accepts as input initial position and velocity vectors (x and y
% coordinates) for a given body. Also accepts M for the central mass in a 
% many-bodied system, assuming M >> m and M >> Nm. Then Sets initial 
% parameters of the two body system and calls ode45, which uses the Runge-
% Kutta method of numerical integration.

% vectors (y0 = [r_x; r_y; v_x; v_y]
y0 = [pos(1), pos(2), vel(1), vel(2)];
time = 0:0.1:time;
error_tol = err;
G = 6.673 * 10^(-11);
M = bigmass; 

% Run the ODE numerical integration. Include options for no bottlenecks.
% note t is unused in second line
options = odeset('RelTol',error_tol,'AbsTol',error_tol);
anonfun = @(t,y)(two_body(t, y, G, M));
[T,Y] = ode45(anonfun, time, y0, options); 

    
end
