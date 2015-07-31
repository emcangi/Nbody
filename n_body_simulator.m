function n_body_simulator()
% Simulates N identical bodies in circular orbits

% Global constants. M = mass of planet, say, Saturn.
G = 6.67384 * 10^(-11);
M = 5.683 *  10^(26);
AU = 149.6e9;  
N = 10000;

% Parameters - radius of orbit set, theta ranges over the unit circle
r = 0.5 * AU;
theta = (2*pi - 0) .* rand(N,1); 
v_max = sqrt(G * M / r);


initial_conditions = zeros(N, 4);            % sun's ICs are all 0
masses = 0.01 * ones(N,1);

% Calculate initial conditions in x, y coordinates and fill array
 for body=1:N
     rx0 = r * cos(theta(body));
     ry0 = r * sin(theta(body));
     vx0 = -v_max * sin(theta(body));
     vy0 =  v_max * cos(theta(body));
     initial_conditions(body, 1) = rx0;
     initial_conditions(body, 2) = ry0;
     initial_conditions(body, 3) = vx0;
     initial_conditions(body, 4) = vy0;
 end

% Set a central, immobile mass, such as a planet Saturn
initial_conditions(1, :) = 0;
masses(1) = M;

[time, data] = main(1, masses, initial_conditions);

% Clean up the data and calculate absolute positions
absolute_positions = zeros(length(time), N);

for body=1:N
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

title('Body distances from orbited mass');
xlabel('Time (days)');
ylabel('Position (AU)');
ylim([0 1]);
set(findall(fig,'-property','FontSize'),'FontSize',16);