%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to simulate N identical bodies in circular orbits around a       %
% massive central object such as a planet.                                %
%                                                                         %
% Author: Eryn Cangi for NSF REU @ CIERA, Northwestern University         %
% Date: 15 June 2015 - present                                            %
% Advisor: Dr. Daniel Abrams                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% GLOBAL CONSTANTS ------------------------------------------------------

G = 6.67384e-11;                                % Nm^2 / kg^2
AU = 149.6e9;                                   % 1 AU in meters

%%% PARAMETERS - CHANGE AS NEEDED -----------------------------------------
%{ 
    r: orbital radius
    M: central mass in kg
    N: number of bodies
    tolerance: absolute and relative
    masses: masses of all bodies
    dt: time intervals to take data, in seconds
    years: total length of time to simulate
    min_sep: minimum separation in radians to ensurea good spread 
%}

r = 0.5 * AU;
M = 5.683e26;%2e30;%        % kg
m = 1e24;
N = 10;
tolerance = 1e-4;
masses = m * ones(N,1);
dt = 24*3600;
years = 1;
min_sep = 0.1;             % For best results do not decrease below 0.05

%%% CALCULATE INITIAL CONDITIONS ------------------------------------------
%{ 
    theta: initial angular positions of bodies in radians
    v: magnitude of tangential velocity at orbit with radius r
    initial_conditions: stores positions and velocities in Cartesian coords
%}

theta = generate_nums(N, min_sep);
% This junk is just for forcing a test of the bodies being too close.
% theta(2) = 0;
% theta(3) = 0.0000222816;
initial_conditions = zeros(N, 4);
v = sqrt(G * M / r);

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
 
%{ 
    for testing
    initial_conditions(2, 3) = -1000;
    initial_conditions(2, 4) = 1000;
    initial_conditions(3, 3) = 1000;
    initial_conditions(3, 4) = -1000;
%}
 
%%% FILE FOR STORING FORCE DATA -------------------------------------------
%{
    fileID = fopen('force_vals.txt','a+');
    fprintf(fileID,'Forces on each body [b1, b2...bn] \r\n Central body is leftmost entry');
%}

%%% SET CENTRAL MASS ------------------------------------------------------
initial_conditions(1, :) = 0;
masses(1) = M;

%%% MAIN SIMULATION CALLER AND DATA OUTPUT --------------------------------
[time, data] = main(years, masses, initial_conditions, tolerance, dt);%, fileID);
%fclose(fileID);

%%% CALCULATE MAGNITUDES OF ORBITAL RADII DURING SIMULATION ---------------
norm_orbit_radii = zeros(length(time), N);

% Get x, y positions and their norm
for body=1:N
    x = data(:, 1, body);
    y = data(:, 2, body);
    norm_orbit_radii(:, body) = sqrt(x .^2 + y .^2);
end 

%%% PLOT RADIUS FROM CENTRAL BODY VS. TIME --------------------------------
fig = figure;
subplot(2,1,1);
hold on;

% Plot all bodies with randomly assigned colors
for idx=1:N
    color = [rand(1) rand(1) rand(1)];
    plot(time, norm_orbit_radii(:,idx),'Color',color); 
end

title(sprintf('Orbital radii vs. time (N=%d, M=%.2e)',N,M));
xlabel('Time (days)');
ylabel('Position (AU)');
ylim([0.49 0.51]);
xlim([0,years*365]);

%%% PLOT INITIAL CONDITIONS -----------------------------------------------

subplot(2,1,2)
hold on;

for idx=1:N
    scatter(initial_conditions(idx,1)/AU, initial_conditions(idx,2)/AU); 
end

title(sprintf('Initial Conditions Visualized'));
xlabel('x (AU)');
ylabel('y (AU)');
plot_limit = r/AU + 0.2;
xlim([-plot_limit,plot_limit]);
ylim([-plot_limit,plot_limit]);
set(findall(fig,'-property','FontSize'),'FontSize',14);
%whitebg(fig,'k');      % toggle for black background

