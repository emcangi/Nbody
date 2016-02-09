%{ A simple script to collide two bodies to test the collision routine.%}

G = 6.67384e-11;                                % Nm^2 / kg^2
AU = 149.6e9;                                   % 1 AU in meters

m = 1e24;
N = 2;
tolerance = 1e-2;
masses = m * ones(N,1);
dt = 24*3600;
years = 1;                  % Must be an integer.
min_sep = 0.05;              % For best results do not decrease below 0.05

theta = generate_nums(N, min_sep);

initial_conditions = zeros(N, 4);

rx1_0 = 0;
ry1_0 = 0;
vx1_0 = 4000;
vy1_0 =  0;

rx2_0 = 1000000;
ry2_0 = 0;
vx2_0 = -4000;
vy2_0 =  0;

initial_conditions(1, 1) = rx1_0;
initial_conditions(1, 2) = ry1_0;
initial_conditions(1, 3) = vx1_0;
initial_conditions(1, 4) = vy1_0;

initial_conditions(2, 1) = rx2_0;
initial_conditions(2, 2) = ry2_0;
initial_conditions(2, 3) = vx2_0;
initial_conditions(2, 4) = vy2_0;


%%% MAIN SIMULATION CALLER AND DATA OUTPUT --------------------------------
[time, data] = main(years, masses, initial_conditions, tolerance, dt);

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

title(sprintf('Orbital radii vs. time'));
xlabel('Time (days)');
ylabel('Position (AU)');
ylim([0 1]);
xlim([0,years*365]);

%%% PLOT INITIAL CONDITIONS -----------------------------------------------

subplot(2,1,2)
hold on;
axis equal;

for idx=1:N
    scatter(initial_conditions(idx,1)/AU, initial_conditions(idx,2)/AU); 
end

title(sprintf('Initial Conditions Visualized'));
xlabel('x (AU)');
ylabel('y (AU)');
%plot_limit = r/AU + 0.2;
% xlim([-plot_limit,plot_limit]);
% ylim([-plot_limit,plot_limit]);
set(findall(fig,'-property','FontSize'),'FontSize',14);
%whitebg(fig,'k');      % toggle for black background

