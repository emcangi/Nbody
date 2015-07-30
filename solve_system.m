function [T, raw_data, AUdata] = solve_system(body_ics, timestep, ...
                                                t_steps, masses, N)
% Solves the second order differential equation d^2r/dt^2 = -(GM/r^3)r.
% Accepts the initial conditions for position and velocity and masses as
% column vectors.

% AU given in meters for unit agreement
AU = 149.6e9;   
SCALE = 1 / AU;

% Parameters for time
dt = timestep;
tf = timestep * t_steps;
time = 0:dt:tf;

% Initial conditons of bodies
y0 = body_ics;

% Establish solver and integration method (Runge-Kutta, order 4)
options = odeset('RelTol',1e-4,'AbsTol',1e-4);
anonfunc = @(t,y)(diffeqs(t, y, masses, N));
[T,Y] = ode45(anonfunc, time, y0, options); 

raw_data = Y;
AUdata = Y * SCALE;

