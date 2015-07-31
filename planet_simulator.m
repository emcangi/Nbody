function planet_simulator()
% Specifically designed to simulate the solar system in 2D.

AU = 149.6e9;  

bodies = 9;

% Options for parameter setting, with the planets. For testing. This will
% be replaced later with a bunch of randomly generated info for identical
% or similar objects.
names = ['Sun    '; 'Mercury'; 'Venus  '; 'Earth  '; 'Mars   '; ...
         'Jupiter'; 'Saturn '; 'Uranus '; 'Neptune'];
namelist = cellstr(names);
axis_set = AU .* [0, 0.387098, 0.723331, 1.00000011, 1.523662, ...
            5.203363, 9.537070, 19.191263, 30.068963];
e_set = [0, 0.20563069, 0.00677323, 0.01671022, 0.09341233, ...
         0.04839266, 0.0541506, 0.04716771, 0.00858587]; 
mass_set = 10^(24) .* [2e6; 0.330; 4.87; 5.97; 0.642; 1898; 568; 86.8; 102];
velocities = 10^3 .* [0, 47.4, 35.0, 29.8, 24.1, 13.1, 9.7, 6.8, 5.4];

initial_conditions = zeros(bodies, 4);            % sun's ICs are all 0
masses = mass_set(1:bodies);

 for planet=2:bodies
     a = axis_set(planet);                    % semi-major axis, AU
     e = e_set(planet);                       % eccentricity
     rx0 = a * (1 + e);                       % Planet initial position
     vy0 = velocities(planet);                % Planet initial velocity
     initial_conditions(planet, 1) = rx0;
     initial_conditions(planet, 4) = vy0;
 end

[time, data] = main(3, masses, initial_conditions);

% Clean up the data and calculate absolute positions
absolute_positions = zeros(length(time), bodies);

for body=1:bodies
    x = data(:, 1, body);
    y = data(:, 2, body);
    absolute_positions(:, body) = sqrt(x .^2 + y .^2);
end 

% Plot all the bodies on one chart
fig = figure(1);
whitebg(1,'k');
hold on;

colors = ['y', 'm', 'y', 'b', 'r', 'r', 'y' 'c', 'b'];

for i=1:9
    plot(time, absolute_positions(:,i), colors(i));
end

title('Planet distances from sun over three years');
xlabel('Time (days)');
ylabel('Position (AU)');
legend(namelist(1:bodies));
set(findall(fig,'-property','FontSize'),'FontSize',16);