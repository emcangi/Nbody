function dy = diffeqs(~, y, masses, N)

% Returns the next state in the differential equation governing motion.
% Elements are: v_x, v_y, GM/r^3 * ?x, GM/r^3 * ?y. MAY BE WRONG.
% y format: [ [rx1, ry1, vx1, vy1];
%             ...[rxN, ryN, vxN, vyN] ] <-- ACTUALLY ITS FLATTENED

% Create column vectors of dr/dt terms where each row is one of the masses
drx_dt = y((N*2)+1 : (N*2)+N);
dry_dt = y((N*3)+1 : (N*3)+N);

% since ICs are passed in as vector, reshape to use this custom function.
% Reshapes into a 
y_reshaped = reshape(y, [N, 4]);
[dvx_dt, dvy_dt] = rhs_of_dvdt(y_reshaped, masses, N);

% put all the data together in one matrix where each row represents data
% for one object, then flatten it back
data_matrix = [drx_dt, dry_dt, dvx_dt, dvy_dt];
flat_result = reshape(data_matrix, [1, N * 4]);

dy = flat_result'; 