clc; clear; close all;
% center of hexagon
x_c = 0;
y_c = 0;

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
figure('Name', '1-1');
hold on;
scatter(0, 0, 100, 'red', 's', 'filled');  % plot the central base station
scatter(x_rand, y_rand, 50, 'blue', 'r', 'filled');

% plot hexagon
plot([x x(1)], [y y(1)], 'k');

xlabel('Distance (m)');
ylabel('Distance (m)');
legend('Central Base Station', 'Mobile Devices', 'Cell Region');
title('Location of Central Base Station and Mobile Devices');
axis square;
grid on;
hold off;



% 六边形半径
L = 1;

% 绘制19个六边形
for i = 1:19
    % 计算六边形的中心点坐标
    if i <= 4
        x_c = i;
        y_c = 1.5;
    elseif i <= 9
        x_c = i - 4.5;
        y_c = 3;
    elseif i <= 15
        x_c = i - 9;
        y_c = 4.5;
    else
        x_c = i - 14.5;
        y_c = 6;
    end
    % 计算六边形的顶点坐标
    x = x_c + L*cosd(0:60:360);
    y = y_c + L*sind(0:60:360);
    % 绘制六边形
    plot([x x(1)], [y y(1)], 'k');
    % 在六边形内随机生成50个点
    x_rand = zeros(50, 1);
    y_rand = zeros(50, 1);
    count = 0;
    while count < 50
        x_rand_temp = rand * 2 * L - L + x_c;
        y_rand_temp = rand * 2 * L - L + y_c;
        if inpolygon(x_rand_temp, y_rand_temp, x, y)
            count = count + 1;
            x_rand(count) = x_rand_temp;
            y_rand(count) = y_rand_temp;
        end
    end
    % 绘制随机点
    hold on;
    scatter(x_rand, y_rand, 10, 'blue', 'filled');
end

% 设置图形属性
axis equal;
axis off;
