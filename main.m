function [time, body_data] = main(years, mass, ic_array, tol)
% Main logic for solving the N-body problem. Receives initial conditions
% for bodies, number of years to simulate. Outputs a tidy 3D array of
% results for each body over time. The format of the results is:
% rows x columns x sheets = time steps x body parameters x bodies.

% Initial Conditions
ics = ic_array;                         % Column vector of row vectors
masses = mass;                          % Masses of bodies: Column vector
N = size(ics, 1);                       % Number of bodies
dt = 24*3600;                           % Time step = 1 earth day
days = 365 * years;                  % Set to 1 less than a year's num
                                        % of days due to indexing troubles

% Flatten the initial conditions array to get 1D column vector.
% format: [x1, x2,...xn, y1...yn, vx1...vxn, vy1...vyn]^T
ics = reshape(ics, [1, N * 4])';

% Do the actual work. results are in the following format:
% [x1, x2...xn, y1, y2...yn, vx1, vx2...vxn, vy1, vy2...vyn]
[T, ~, AUdata] = solve_system(ics, dt, days, masses, N, tol);

% Makes a matrix where: rows = times, columns = values, pages = bodies
rows = size(T,1);
AUdata_tidy = zeros(rows, 4, N);

for m=1:rows
    temp = reshape(AUdata(m, :), [N,4]);
    for n=1:N
        AUdata_tidy(m, :, n) = temp(n, :);
    end
end

% Results
time = T * (1/dt);                       % Convert T array to days
body_data = AUdata_tidy;


