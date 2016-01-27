function [dvx_dt, dvy_dt] = rhs_of_dvdt(y, masses, N)%, fileID)
%{ 
    Calculates GMd / d^3 component-wise between each mass and each other
    mass
    y: a flat array of initial conditions, form rx1, ry1, vx1, vy1, rx2..
    masses: N x 1 array of mass values for each body
    N: number of bodies
%}

G = 6.67384 * 10^(-11);

%%% GET POSITIONS AND CALCULATE ALL SEPARATIONS ---------------------------

% Make new column vectors with only the positions (x and y each)
x_positions = y(:, 1);
y_positions = y(:, 2);

%{ 
    Create an "each to each" set of differences between objects
    format is like: [x1-x1, x1-x2, ... x1-xN; x2-x1...x2-xN; xN-x1...xN-xN]
%}
x_diffs = repmat(x_positions', N, 1) - repmat(x_positions, 1, N);
y_diffs = repmat(y_positions', N, 1) - repmat(y_positions, 1, N);

%%% CALCULATE GMd/d^3 ------------------------------------------------------

%{ 
    Replicate the masses column vector N times to use in element-wise
    calculation
%}
m = repmat(masses, 1, N);

%{
    calculate G * other mass * distance along an axis, ELEMENT WISE.
    Therefore, each column corresponds to the calculations for one body.
%}
GMrx = -G * (m .* x_diffs);
GMry = -G * (m .* y_diffs);

% Square differences element-wise
x_sqrs = x_diffs .^ 2;
y_sqrs = y_diffs .^ 2;

%{ 
    Take the sum of the squares, then take square root element-wise, then
    cube the answer. This gives an array that holds the cube of the mag of
    the distances between each two bodies, like:
    [0, d2-1, d3-1; d1-2, 0, d3-2; ...]
%}
sum_sqrs = x_sqrs + y_sqrs;
r = sqrt(sum_sqrs);
d3 = r .^ 3;

%{ 
    Divide by d^3, ELEMENT WISE. In this result, each column corresponds to
    all the GMd_x/d^3 and GMd_y/d^3 terms for one object. So the first
    column is the column where m1 is the "self" mass and so on.
%}
all_dvx_dt = GMrx ./ d3;
all_dvy_dt = GMry ./ d3;

% Handle NaN caused by dividing by 0 by 0 when self and other are the same
all_dvx_dt(isnan(all_dvx_dt)) = 0;
all_dvy_dt(isnan(all_dvy_dt)) = 0;

% Sum each column to get the collective acceleration on each mass
net_accel_x = sum(all_dvx_dt);
net_accel_y = sum(all_dvy_dt);

%%% CHECKING THE NET FORCES -----------------------------------------------
%{ 
    This section is for troubleshooting. It is not necessary for the
    simulation to run.

    Because all the masses are currently the same, I can just multiply the
    whole acceleration array by 1e24 to find the forces.

    net_F_x = 1e24 .* net_accel_x;
    net_F_y = 1e24 .* net_accel_y; 
    net_F = sqrt(net_F_x .^ 2 + net_F_y .^ 2);
    string = mat2str(net_F);
    fprintf(fileID, '%s\r\n',string);
%}

%%% RESULTS dv/dt ---------------------------------------------------------

%{ 
    Transpose to get column vectors where each row has a single value for
    each object: [body 1 accel; body 2 accel; ... body N accel]
%}
dvx_dt = net_accel_x.';
dvy_dt = net_accel_y.';

end







