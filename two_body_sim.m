AU = 149.6e9;
SCALE = 250 / AU;
angle = sqrt(2)/2;

% simulate earth and sun
pos = [AU*angle, AU*angle];
vel = [3e4*angle, 3e4*angle];
err = 1e-3;
bigmass = 2e30;
time = 10000;

[T, Y] = two_body_caller(pos, vel, err, bigmass, time);