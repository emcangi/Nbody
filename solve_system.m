function [T_data, raw_data, AUdata] = solve_system(body_ics, timestep, ...
                                              days, masses, N, tol)
% Solves the second order differential equation d^2r/dt^2 = -(GM/r^3)r
% Where r is a vector. Accepts the initial conditions for position and 
% velocity and masses as column vectors.

% Variables and parameters
global dists;
AU = 149.6e9;   
SCALE = 1 / AU;
ti = 0;
dt = timestep;
tf = dt * (days - 1);                % -1 makes row num right in time array
y0 = body_ics;
t_out = zeros(days+1, 1);            % hold accumulated times
y_out = zeros(days+1, N * 4);        % hold accumulated data
y_out(1, :) = y0';                      
fill_start = 2;                     % start filling after first entry 
                                    % (which is initial conditions)
in_progress = true;

% Establish solver and integration method (Runge-Kutta, order 4)
options = odeset('RelTol', tol, 'AbsTol', tol, 'Events', @events);
anonfunc = @(t,y)(diffeqs(t, y, masses, N));

% Runs the integration. Events cause it to stop, so it is in a loop. 
while in_progress
    time = ti:dt:tf;
    [T,Y,~,YE,IE] = ode45(anonfunc, time, y0, options); 

    if ~isempty(YE)
        % runs if YE has values, i.e. an event occurred
        disp('Collision!');
        nt = length(T) - 1;             % successful steps minus ICs
        t_out(fill_start:fill_start + nt -1) = T(2:nt +1);    
        y_out(fill_start:fill_start + nt -1, :) = Y(2:nt+1, :);
        event_value = dists(IE);

        % Dists was originally in an upper triangular matrix. This
        % block converts it back with sneaky tricks.
        temp = tril(ones(N),-1);
        temp(~~temp) = dists;
        dists_ut = temp';

        % Find the row and column of the event value. These each
        % represent a body that was involved in the collision.
        [b1, b2] = find(dists_ut==event_value);

        % Collision handling and velocity reassignment
        vxi1 = YE(3 + 4*(b1-1));          % body1 x velocity, initial
        vyi1 = YE(4 + 4*(b1-1));          % body1 y velocity, initial
        vxi2 = YE(3 + 4*(b2-1));          % body2 x velocity, initial
        vyi2 = YE(4 + 4*(b2-1));          % body2 y velocity, initial

        [vxf1, vyf1, vxf2, vyf2, fail] = collision(vxi1, vyi1, vxi2, vyi2); % REPLACE
      
        if ~fail
            disp('Successful collision!')
            YE(3 + 4*(b1-1)) = vxf1;          % body1 x velocity, after
            YE(4 + 4*(b1-1)) = vyf1;          % body1 y velocity, after
            YE(3 + 4*(b2-1)) = vxf2;          % body2 x velocity, after
            YE(4 + 4*(b2-1)) = vyf2;          % body2 y velocity, after
            y0 = YE;                            % reset initial conditions
        end
        
        % NEW TIME START HERE - might need to fix time step
        ti = T(nt);
        fprintf('New start time: day %f \n', ti/timestep)
        fill_start = fill_start + nt;
    end

    % Finished integrating if most recent last time step is equal to
    % the designated ending time step
    if T(end) >= tf
        nt = length(T) - 1;             % successful steps minus ICs
        t_out(fill_start:end) = T(2:nt);    
        y_out(fill_start:end, :) = Y(2:nt, :);
        break
    end
end

% Determine which array has the final data. 
if length(T) >= days                % No event occurred; T is cumulative
    raw_data = Y;
    T_data = T;
else                                   % Event occurred
    raw_data = y_out;
    T_data = t_out;
end

AUdata = raw_data * SCALE;

function [value, isterminal, direction] = events(~,y)
    % Event function to check if distances between any two particles
    % are under a certain threshhold (1000 km).     

    rx = y(1:N);
    ry = y(N+1:2*N);
    N_separations = N * (N-1) / 2;

    % Find all distances between particles on each axis. 3rd and 6th
    % lines use sneaky tricks to flatten upper triangular matrix into a
    % column vector in which elements are filled by reading the upper
    % triangular matrix left to right, top to bottom. 
    [x_pos1, x_pos2] = meshgrid(rx);
    x_diffs = x_pos1 - x_pos2;
    x_dists = x_diffs(tril(true(size(x_diffs)), -1));
    [y_pos1, y_pos2] = meshgrid(ry);
    y_diffs = y_pos1 - y_pos2;
    y_dists = y_diffs(tril(true(size(y_diffs)), -1));

    % Calculate distances and convert to column vector with nonzeros()
    dists = sqrt(x_dists .^ 2 + y_dists .^ 2);

    % detect when objects are within 1 km of each other and stop if so
    % MAY NOT WORK IF THERE ARE MORE THAN ONE COLLISION? CHECK LATER
    value = dists - 100000;
    isterminal = ones(N_separations, 1);
    direction = -1 * ones(N_separations, 1);
end
end
