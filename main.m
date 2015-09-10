function [time, body_data] = main(years, masses, ic_array, tol, dt)
% Main logic for solving the N-body problem. Receives: masses, initial
% conditions for all bodies, number of years to simulate and a tolerance
% value. Outputs a 3D array of position and velocity data for each body
% over time. Rows represent steps in time. Columns represent data values in
% this format: x position, y position, x velocity, y velocity. Each sheet
% represents one body in the system.

% INITIAL CONDITIONS ------------------------------------------------------
N = size(ic_array, 1);
days = 365 * years;

% Solver only accepts flat arrays, so this line flattens the initial
% conditions array to get 1D column vector.
% format: [x1; x2;...xn; y1...yn; vx1...vxn; vy1...vyn]
ics = reshape(ic_array, [1, N * 4])';

% SOLVE THE SYSTEM --------------------------------------------------------
% Result format: [x1; x2;...xn; y1...yn; vx1...vxn; vy1...vyn]
% To return raw data in kg/m/s units, rename the ~ to whatever you want
[T, ~, AUdata] = solve_system(ics, dt, days, masses, N, tol);


% UNFLATTEN THE OUTPUT DATA, MAKE OUTPUT ARRAY ----------------------------
rows = size(T,1);
AUdata_3d = zeros(rows, 4, N);

for m=1:rows
    temp = reshape(AUdata(m, :), [N,4]);
    for p=1:N
        AUdata_3d(m, :, p) = temp(p, :);
    end
end

% RESULTS -----------------------------------------------------------------
time = T * (1/dt);                                % Convert T array to days
body_data = AUdata_3d;

end
