function dy = diffeqs(~, y, masses, N)%, fileID)
%{ 
    Solves the second order differential equation d^2r/dt^2 = -(GM/d^3)*d
    where r is a positional vector and d is the separartion between any two
    objects. Uses coupled equations dr/dt = v and dv/dt = -GM/d^3 * r where
    v and r are vectors.
    y format (Required): [rx1...rxN; ry1...ryN; vx1...vxN; vy1...vyN]
%}

%%% dr/dt = v -------------------------------------------------------------
drx_dt = y((N*2)+1 : (N*2)+N);
dry_dt = y((N*3)+1 : (N*3)+N);

%%% dv/dt = -GM/d^3 * d ---------------------------------------------------
%{ 
    Reshapes into format: [[rx1, ry1, vx1, vy1]; ...[rxN, ryN, vxN, vyN]]
    which is a bit easier to understand and used within rhs_of_dvdt
    function 
%}

y_reshaped = reshape(y, [N, 4]);
[dvx_dt, dvy_dt] = rhs_of_dvdt(y_reshaped, masses, N);%, fileID);

%{ 
    put all the data together in data_matrix where each row holds x, y, vx
    and vy for one object, then flatten it again
%}
data_matrix = [drx_dt, dry_dt, dvx_dt, dvy_dt];
flat_result = reshape(data_matrix, [1, N * 4]);

dy = flat_result'; 

end