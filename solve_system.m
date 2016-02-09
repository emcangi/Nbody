function [T_data, raw_data, AUdata]=solve_system(y0, dt, days, mass, N, ...
                                                 tol)%, fileID)
%{ 
    Solves the differential equation, handles collision events and places
    output data into arrays. 
    Accepts as Input:
    y0: initial conditions
    dt:  points in time at which to record data
    days: number of days to run simulation
    mass: masses of all objects
    N: number of bodies
    tol: error tolerance, both relative and absolute
%}

%%% VARIABLES AND PARAMETERS ----------------------------------------------
%{
    dists contains separations between each object and each other object.
    It must remain a global variable so it can be acquired from the event 
    handler; it is used to identify which objects had collisions.
%}
global dists;

AU = 149.6e9;   
SCALE = 1 / AU;
ti = 0;
tf = dt * (days - 1);               % -1 ensures correct number of rows in 
                                    % time array
  
%{
   t_out and y_out serve the same purpose as the solver outputs T and Y,
   but these are used in simulations where collisions occur to manually
   keep track of the data since events stop the simulation and require
   restarting                                 
%}
t_out = zeros(days, 1);             % add 1 because ti = 0 not ti = 1
y_out = zeros(days, N * 4);
y_out(1, :) = y0';
marker = 1;                         % used to track where data left off 
                                    % after an event. first row is ICs
in_progress = true;

%%% SET SOLVER AND INTEGRATION METHOD -------------------------------------
% Using an anonymous function allows passing parameters to diffeqs
options = odeset('RelTol', tol, 'AbsTol', tol, 'Events', @events);
anonfunc = @(t,y)(diffeqs(t, y, mass, N));%, fileID));

%%% RUN THE SIMULATION ----------------------------------------------------
% Looped due to necessity of stopping and restarting after events
while in_progress
    time = ti:dt:tf;
    [T,Y,~,YE,IE] = ode45(anonfunc, time, y0, options); 

    % runs if YE has values, i.e. an event (collision) occurred
    if ~isempty(YE)
        s = length(T) - 1;          % Saveable data, excluding collisions
        ti = T(s) + 86400;          % Restart index day after last complete
                                    % before most recent collision event
        currday = T(s)/dt;          % number of current day
        alert = sprintf('Collision on day %d',currday);
        alert2 = sprintf('Index of collision in meshgrid: %d',IE);
        disp(alert);
        disp(alert2);
        t_out(marker:marker + s - 1) = T(1:s);
        y_out(marker:marker + s - 1, :) = Y(1:s, :);
        
        event_value = dists(IE);    % find which bodies had collision

        %{ 
            In event handler, dists is an upper triangular matrix. This
            block converts it back so we can find the right bodies. See
            event handler for details on the format of this.
        %}
        unflat_YE = tril(ones(N),-1);
        unflat_YE(~~unflat_YE) = dists;
        dists_ut = unflat_YE';

        %{
            Find the row and column of the event value. These each
            represent a body that was involved in the collision.
        %}
        [b1, b2] = find(dists_ut==event_value);
        
        fprintf('Colliding bodies: %d, %d',b1,b2);
        fprintf('\n');

        % Create an unflattened copy of YE so it's easier to grab the right
        % bodies
        YE_unflat = reshape(YE(1, :), [N,4]);
        
        % Collect pre-collision data for one of the bodies
        rxi1 = YE_unflat(b1, 1);
        ryi1 = YE_unflat(b1, 2);
        vxi1 = YE_unflat(b1, 3);
        vyi1 = YE_unflat(b1, 4);
        
        % Collect pre-collision data for the other of the bodies
        rxi2 = YE_unflat(b2, 1);
        ryi2 = YE_unflat(b2, 2);
        vxi2 = YE_unflat(b2, 3);
        vyi2 = YE_unflat(b2, 4);
        
        % Print stuff for testing
        fprintf('Precollision data:\n');
        fprintf('Body %d position: (%.2e, %.2e)\n',b1,rxi1,ryi1);
        fprintf('Body %d velocity: (%.2e, %.2e)\n',b1,vxi1,vyi1);
        fprintf('Body %d position: (%.2e, %.2e)\n',b2,rxi2,ryi2);
        fprintf('Body %d velocity: (%.2e, %.2e)\n',b2,vxi2,vyi2);
        
        % Calculate post-collision velocities and a success flag
        [vxf1, vyf1, vxf2, vyf2, fail] = collision(rxi1, ryi1, rxi2, ryi2, ...
                                        vxi1, vyi1, vxi2, vyi2, dt);
                                    
        fprintf('Postcollision data:\n');
        fprintf('Body %d velocity: (%.2e, %.2e)\n',b1,vxf1,vyf1);
        fprintf('Body %d velocity: (%.2e, %.2e)\n',b2,vxf2,vyf2);
      
        % Assign new velocities
        if ~fail
            disp('Successful collision!')
            disp('')
            YE_unflat(b1,3) = vxf1; % body1 x velocity, after
            YE_unflat(b1,4) = vyf1; % body1 y velocity, after
            YE_unflat(b2,3) = vxf2; % body2 x velocity, after
            YE_unflat(b2,4) = vyf2; % body2 y velocity, after
            
            % Reflatten YE
            y0 = reshape(YE_unflat, [1, N * 4])'; % reset initial conditions
        end
        
        % Update the marker that shows next fill position in array
        fprintf('New start time: day %d \n',currday+1)
        marker = marker + s;
    end

    %{ 
        Find out if finished with all time steps and collect data that
        occurred after the last collision
    %}
    if T(end) >= tf
        t_out(marker:end) = T(1:end);    
        y_out(marker:end, :) = Y(1:end, :);
        break
    end
end

%%% ASSIGN WHICHEVER DATA ARRAY IS COMPLETE TO OUTPUT ARRAY ---------------
if length(T) >= days              % No event occurred and T had max entries
    raw_data = Y;
    T_data = T;
else                              % Event occurred; use the managed data
    raw_data = y_out;
    T_data = t_out;
end

% Scale the data to be in AUs instead of meters
AUdata = raw_data * SCALE;

%%% EVENT HANDLER ---------------------------------------------------------
function [value, isterminal, direction] = events(~,y)
    %{ 
        Event function to check if distances between any two particles are
        under a certain threshhold.     
    %}
    
    threshold = 1e6;
    
    % Gather the position data from the flattened array
    rx = y(1:N);
    ry = y(N+1:2*N);
    pairs = N * (N-1) / 2;      % number of unique object pairs

    %{ 
        Find all distances between objects on each axis. 3rd and 6th
        lines use sneaky tricks to flatten upper triangular matrix into a
        column vector in which elements are filled by reading the upper
        triangular matrix left to right, top to bottom. 
    %}
    [x_pos1, x_pos2] = meshgrid(rx);
    x_diffs = x_pos1 - x_pos2;
    x_dists = x_diffs(tril(true(size(x_diffs)), -1));
    
    [y_pos1, y_pos2] = meshgrid(ry);
    y_diffs = y_pos1 - y_pos2;
    y_dists = y_diffs(tril(true(size(y_diffs)), -1));

    %{
        Calculate distances and convert to column vector, required for
        logic of event handling below
    %}
    dists = sqrt(x_dists .^ 2 + y_dists .^ 2);

    % *** MAY NOT WORK IF THERE ARE MORE THAN ONE COLLISION? CHECK LATER
    value = dists - threshold;
    isterminal = ones(pairs, 1);
    direction = -1 * ones(pairs, 1);
end
end
