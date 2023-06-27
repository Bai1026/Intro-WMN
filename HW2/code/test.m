clc;clear;close all;

num_points = 50;
L = 1;
x_c = 0;
y_c = 0;

% plot the first hexagon
figure('Name', 'Cellular Network');
hold on;

for i = 1:19
    % hexagon's vertices
    x = x_c + L*cosd(0:60:360);
    y = y_c + L*sind(0:60:360);

    % generate 50 points that's in the hexagon
    x_rand = zeros(num_points, 1);
    y_rand = zeros(num_points, 1);
    count = 0;

    while count < num_points
        x_rand_temp = rand * 2 * L - L;
        y_rand_temp = rand * 2 * L - L;
        if inpolygon(x_rand_temp, y_rand_temp, x, y)
            count = count + 1;
            x_rand(count) = x_rand_temp;
            y_rand(count) = y_rand_temp;
        end
    end

    
    % plot the points in the hexagon
    scatter(x_rand, y_rand, 50, 'blue', 'r', 'filled');

    % plot hexagon
    plot([x x(1)], [y y(1)], 'k');

    % update center point for next iteration
    x_c = x_c + 1.5;
    y_c = y_c - sqrt(3)/2 + 2*sqrt(3);
end

% plot the central base station
scatter(0, 0, 100, 'red', 's', 'filled');

xlabel('Distance (m)');
ylabel('Distance (m)');
legend('Mobile Devices', 'Cell Region', 'Central Base Station');
title('Location of Central Base Station and Mobile Devices');
axis square;
grid on;
hold off;
