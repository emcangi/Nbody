function [T, raw_data, AUdata] = solve_system(body_ics, timestep, ...
                                              t_steps, masses, N)
    % Solves the second order differential equation d^2r/dt^2 = -(GM/r^3)r
    % Where r is a vector. Accepts the initial conditions for position and 
    % velocity and masses as column vectors.

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
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-4, 'Events', @events);
    anonfunc = @(t,y)(diffeqs(t, y, masses, N));
    [T,Y,TE,YE,IE] = ode45(anonfunc, time, y0, options); 
    
    % do something with TE, YE, IE
    disp(TE);
    
    raw_data = Y;
    AUdata = Y * SCALE;

function [value, isterminal, direction] = events(t,y)
    
    rx = y(1:N);
    ry = y(N+1:2*N);
    N_separations = N * (N-1) / 2;
    
    % Use meshgrid again to find all distances between particles on each
    % axis. Then use triu to remove duplicates and nonzeros to flatten into
    % a column vector
    [x_pos1, x_pos2] = meshgrid(rx);
    x_diffs = x_pos1 - x_pos2;
    x_dists = nonzeros(triu(x_diffs, 1));
    [y_pos1, y_pos2] = meshgrid(ry);
    y_diffs = y_pos1 - y_pos2;
    y_dists = nonzeros(triu(y_diffs, 1));   
    
    % Find the square of the distances between particles so we only have
    % one value to work with
    dist_squared = x_dists .^ 2 + y_dists .^ 2;
    
    % detect when objects are within 1 km of each other
    value = dist_squared - 1000^2;
    isterminal = ones(N_separations, 1);
    direction = -1 * ones(N_separations, 1);
end
end
