% Plot charts of positions

[names, time, data] = main(3);

sun_x = data(:, 1, 1); 
sun_y = data(:, 2, 1);
sun_positions = sqrt(sun_x .^2 + sun_y .^2);

earth_x = data(:, 1, 2);
earth_y = data(:, 2, 2);
% inverse = earth_x .^ -1;
% earth_sign = earth_x .* abs(inverse)
earth_positions = sqrt(earth_x.^2 + earth_y.^2);

mars_x = data(:, 1, 3);
mars_y = data(:, 2, 3);
% inverse = mars_x .^ -1;
% mars_sign = mars_x .* abs(inverse)
mars_positions = sqrt(mars_x .^2 + mars_y .^2);

fig = figure(1);
whitebg(1,'k')
plot(time(1:5:end), sun_positions(1:5:end), 'yo', ...
     time(1:5:end), earth_positions(1:5:end), 'bo', ...
     time(1:5:end), mars_positions(1:5:end), 'r.', ...
     time, sun_positions, 'y', ...
     time, earth_positions, 'b', ...
     time, mars_positions, 'r');
title('Distance from origin (sun) for Earth and Mars over three years');
xlabel('Time (days)');
ylabel('Position (AU)');
legend('Sun', 'Earth', 'Mars');
set(findall(fig,'-property','FontSize'),'FontSize',16) 