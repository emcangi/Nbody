function [time, body_data] = main(years, masses, ic_array, tol, dt)%, fileID)
%{
    This script does a little bit of housekeeping on data sent to and
    received from the solve_system.m function, which provides the main
    logic of the simulation. Mostly, this script deals with reshaping
    arrays of positions and velocities since the solver requires a flat
    array but that is really hard to read.

    Receives: simulation length, masses, initial conditions, tolerance,
    timestep and a filename if needed.

    Outputs: 3D array of position and velocity over time. Rows = entries in
    time; columns = data values in format x position, y position, x
    velocity, y velocity; sheet = one body each.
%}  

%%% INITIAL CONDITIONS ----------------------------------------------------
N = size(ic_array, 1);
days = 365 * years;

% Flattens the initial conditions array to get 1D column vector.
% format: [x1; x2;...xn; y1...yn; vx1...vxn; vy1...vyn]
ics = reshape(ic_array, [1, N * 4])';


%%% SOLVE THE SYSTEM ------------------------------------------------------
% Result format is flat: [x1; x2;...xn; y1...yn; vx1...vxn; vy1...vyn]
% To return raw data in kg/m/s units, rename the ~ to whatever you want
[T, ~, AUdata] = solve_system(ics, dt, days, masses, N, tol);%, fileID);


%%% UNFLATTEN THE OUTPUT DATA, MAKE OUTPUT ARRAY --------------------------
rows = size(T,1);
AUdata_3d = zeros(rows, 4, N);

for m=1:rows
    temp = reshape(AUdata(m, :), [N,4]);
    for p=1:N
        AUdata_3d(m, :, p) = temp(p, :);
    end
end

%%% RESULTS ---------------------------------------------------------------
time = T * (1/dt);                                % Convert T array to days
body_data = AUdata_3d;

end
