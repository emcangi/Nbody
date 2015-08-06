% Specifically designed to simulate the solar system in 2D.

AU = 149.6e9;  
G = 6.67384 * 10^(-11);
bodies = 5;

% Options for parameter setting, with the planets. Axes are the semimajor
% axis of each planet, and peri represents perihelion.
names = ['Sun    '; 'Mercury'; 'Venus  '; 'Earth  '; 'Mars   '; ...
         'Jupiter'; 'Saturn '; 'Uranus '; 'Neptune'];
namelist = cellstr(names);
peri_set = AU .* [0, 0.30749951, 0.7184327, 0.98328989, 1.388133346, ...
                  4.95155843, 9.020632, 18.286056, 29.810795];
axis_set = AU .* [0, 0.387098, 0.723331, 1.00000011, 1.523662, ...
                  5.203363, 9.537070, 19.191263, 30.068963];
%e_set = [0, 0.20563069, 0.00677323, 0.01671022, 0.09341233, ...
%         0.04839266, 0.0541506, 0.04716771, 0.00858587]; 
mass_set = 10^(24) .* [2e6; 0.330; 4.87; 5.97; 0.642; 1898; 568; 86.8; 102];

initial_conditions = zeros(bodies, 4);            % sun's ICs are all 0
masses = mass_set(1:bodies);

 for planet=2:bodies
     a = axis_set(planet);
     velocity = sqrt(G*masses(1)*(2/peri_set(planet) - 1/a));
     %e = e_set(planet);                      % eccentricity
     rx0 = peri_set(planet);                  % Planet initial position
     vy0 = velocity;
     initial_conditions(planet, 1) = rx0;
     initial_conditions(planet, 4) = vy0;
 end

[time, data] = main(3, masses, initial_conditions);

% Clean up the data and calculate absolute positions
absolute_positions = zeros(length(data), bodies);

for body=1:bodies
    x = data(:, 1, body);
    y = data(:, 2, body);
    absolute_positions(:, body) = sqrt(x .^2 + y .^2);
end 

% Plot all the bodies on one chart
fig = figure(1);
whitebg(1,'k');
hold on;

colors = ['g', 'm', 'y', 'b', 'r', 'r', 'y' 'c', 'b'];

for i=1:bodies
    plot(time, absolute_positions(:,i), colors(i));
end

title('Planet distances from sun over three years');
xlabel('Time (days)');
ylabel('Position (AU)');
legend(namelist(1:bodies));
set(findall(fig,'-property','FontSize'),'FontSize',16);