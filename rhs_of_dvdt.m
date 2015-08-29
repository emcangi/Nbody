function [dvx_dt, dvy_dt] = rhs_of_dvdt(y, masses, N)
% Calculates GMr_x / r^3 and GMr_y / r^3 between each mass and each other
% mass, using repmat and careful control of element-wise and 
% non-element-wise operations.

G = 6.67384 * 10^(-11);

% Make new column vectors with only the positions (x and y each)
x_positions = y(:, 1);
y_positions = y(:, 2);

% Create an "each to each" set of differences among particle locations
x_diffs = repmat(x_positions', N, 1) - repmat(x_positions, 1, N);
y_diffs = repmat(y_positions', N, 1) - repmat(y_positions, 1, N);

% Replicate the masses column vector N times to use in element-wise
% calculation
m = repmat(masses, 1, N);

% calculate G * other mass * distance along an axis, ELEMENT WISE. The
% result is like:
% [ 0, m1(rx2 - rx1), m1(rx3 - rx1);
%   m2(rx1 - rx2), 0, m2(rx3 - rx2); ...]
% Therefore, each column corresponds to the calculations for one body.
GMrx = -G * (m .* x_diffs);
GMry = -G * (m .* y_diffs);

% Now we can square the differences element-wise
x_sqrs = x_diffs .^ 2;
y_sqrs = y_diffs .^ 2;

% Take the sum of the squares, then take square root element-wise, then
% cube the answer. This gives an array that holds the cube of the mag of
% the distances between each two bodies, like:
% [0, r2-1, r3-1; r1-2, 0, r3-2; ...]
sum_sqrs = x_sqrs + y_sqrs;
r = sqrt(sum_sqrs);
r3 = r .^ 3;

% Divide by r^3, ELEMENT WISE. In this result, each column corresponds to
% all the GMrx/r^3 and GMry/r^3 terms for a mass. So the first column is
% the column where m1 is the "self" mass and so on.
all_dvx_dt = GMrx ./ r3;
all_dvy_dt = GMry ./ r3;

% handle NaN caused by dividing by 0 by 0
all_dvx_dt(isnan(all_dvx_dt)) = 0;
all_dvy_dt(isnan(all_dvy_dt)) = 0;

% Now we can sum each column to get the collective acceleration on each
% mass, and finally, transpose again so we get back column vectors.
net_accel_x = sum(all_dvx_dt);
net_accel_y = sum(all_dvy_dt);

% these are the column vectors where each row contains the net acceleration
% on a single body.
dvx_dt = net_accel_x.';
dvy_dt = net_accel_y.';







