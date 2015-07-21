function dy = two_body(~, y, G, M)
% Accepts as input y = [r; v] which are the initial conditions and can be 
% expanded to components as [r_x, r_y, v_x, v_y]. Also accepts G and M
% which are constants, G = Gravitational constant, M = central mass

dy = (1:4)';

% dy = [drx/dt, dry/dt, dvx/dt, dvy/dt]
dy(1) = y(3);                           % dy(1) = v_x
dy(2) = y(4);                           % dy(2) = v_y
dy(3) = -(G*M)/((y(1))^3);                % dy(3) = -(G*M)/(rx^3)
dy(4) = -(G*M)/((y(2))^3);                % dy(4) = -(G*M)/(ry^3)

end