%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to simulate N identical bodies in circular orbits around a       %
% massive central object such as a planet.                                %
%                                                                         %
% Author: Eryn Cangi for NSF REU @ CIERA, Northwestern University         %
% Date: 15 June 2015 - present                                            %
% Advisor: Dr. Daniel Abrams                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GLOBAL CONSTANTS --------------------------------------------------------
% M = mass of planet, e.g. Saturn. AU = in meters, N = number of bodies. 
% tolerance is applied to both absolute and relative.
G = 6.67384e-11;
M = 2e30;%5.683e26;
AU = 149.6e9;
N = 40;         
tolerance = 1e-2;
hello

% PARAMETERS - CHANGE AS NEEDED -------------------------------------------
% r = orbital radius, theta = locations, v = magnitude of tangential
% velocity at this orbit, initial_conditions stores positions and
% velocities, masses designates orbiting body masses, and dt is the time
% step to examine in seconds. (24*3600 = 1 earth day)

% Random seed for body locations
rng(1);

r = 0.5 * AU;
theta = (2*pi - 0) .* rand(N,1); 
v = sqrt(G * M / r);
initial_conditions = zeros(N, 4);
masses = 1e24 * ones(N,1);
timestep = 24*3600;
years = 3;

% CALCULATE INITIAL CONDITIONS --------------------------------------------
% using trigonometry and Cartesian coordinates
 for body=1:N
     rx0 = r * cos(theta(body));
     ry0 = r * sin(theta(body));
     vx0 = -v * sin(theta(body));
     vy0 =  v * cos(theta(body));
     initial_conditions(body, 1) = rx0;
     initial_conditions(body, 2) = ry0;
     initial_conditions(body, 3) = vx0;
     initial_conditions(body, 4) = vy0;
 end

% SET ORBITED MASS --------------------------------------------------------
initial_conditions(1, :) = 0;
masses(1) = M;

% MAIN SIMULATION CALLER AND DATA OUTPUT ----------------------------------
[time, data] = main(years, masses, initial_conditions, tolerance, timestep);

% CALCULATE MAGNITUDES OF ORBITAL RADII DURING SIMULATION -----------------
norm_orbit_radii = zeros(length(time), N);

% Get x, y positions and their norm
for body=1:N
    x = data(:, 1, body);
    y = data(:, 2, body);
    norm_orbit_radii(:, body) = sqrt(x .^2 + y .^2);
end 

% PLOT RADIUS FROM CENTRAL BODY VS. TIME ----------------------------------
fig = figure;
whitebg(fig,'k');
hold on;

% Plot all bodies with randomly assigned colors
for idx=2:N
    color = [rand(1) rand(1) rand(1)];
    plot(time, norm_orbit_radii(:,idx),'Color',color); 
end

title(sprintf('Body distances from orbited mass, N = %d',N));
xlabel('Time (days)');
ylabel('Position (AU)');
ylim([0 1]);
set(findall(fig,'-property','FontSize'),'FontSize',14);