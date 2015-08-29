% Simulates N identical bodies in circular orbits

rng(2);                                     % random seed

% Global constants. M = mass of planet, say, Saturn.
G = 6.67384e-11;
M = 5.683e26;
AU = 149.6e9;  
N = 40;
tolerance = 1e-5;

% Parameters - radius of orbit set, theta ranges over the unit circle
r = 0.5 * AU;
theta = (2*pi - 0) .* rand(N,1); 
v = sqrt(G * M / r);
initial_conditions = zeros(N, 4);            % sun's ICs are all 0
masses = 1e24 * ones(N,1);

% Calculate initial conditions in x, y coordinates and fill array
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

% Set a central, immobile mass, such as a planet ~Saturn
initial_conditions(1, :) = 0;
masses(1) = M;

% Do the math
[time, data] = main(3, masses, initial_conditions, tolerance);

% Clean up the data and calculate absolute positions
absolute_positions = zeros(length(time), N);

for body=1:N
    x = data(:, 1, body);
    y = data(:, 2, body);
    absolute_positions(:, body) = sqrt(x .^2 + y .^2);
end 

% Plot all the bodies on one chart
fig = figure;
whitebg(fig,'k');
hold on;

% Plot all bodies with randomly assigned colors
for idx=2:N
    color = [rand(1) rand(1) rand(1)];
    plot(time, absolute_positions(:,idx),'Color',color); 
end


title(sprintf('Body distances from orbited mass, N = %d',N));
xlabel('Time (days)');
ylabel('Position (AU)');
ylim([0 1]);
set(findall(fig,'-property','FontSize'),'FontSize',14);