clc; clear; close all;
%------------------------1.1-------------------------
num_cell = 19;
isd = 500;
radius = isd/sqrt(3);
distances = zeros(num_cell);
angles = zeros(num_cell);

index1 = 2:7;
index2 = 8:13;
index3 = 14:19;
distances(index1) = isd;
distances(index2) = 2*isd*cosd(30);
distances(index3) = 2*isd;
angles(index1) = 30:60:360;
angles(index2) = 0:60:300;
angles(index3) = 30:60:360;

% calculate the coordinates of x and y
x_center = 250;
y_center = 0;

all_x = zeros(num_cell*7, 1);
all_y = zeros(num_cell*7, 1);

for i = 1:num_cell
    all_x(i) = x_center + distances(i)*cosd(angles(i));
    all_y(i) = y_center + distances(i)*sind(angles(i));
end

% plot the hexagonal cells
figure;
hold on;
plot_circles_and_labels(radius, all_x, all_y)
plot_all(radius, isd, all_x, all_y)
hold off;

%------------------------1.2-------------------------
% initialize mobile device location, velocity, direction and time interval
x = 250;
y = 0;
velocity = 1 + 14*rand(); % randomly generate velocity between 1 and 15 m/s
direction = 2*pi*rand(); % randomly generate direction between 0 and 2*pi
t_interval = 1 + 5*rand(); % randomly generate time interval between 1 and 6 seconds
t = 0; % initialize time 
handoff_record = []; % initialize handoff record
total_time = 900; % total simulation time

% Create plot of the hexagonal cells
figure;
hold on;
plot_circles_and_labels(radius, all_x, all_y)
plot_all(radius, isd, all_x, all_y)

% plot the initial position of mobile device
plot(x, y, 'o-', 'MarkerSize', 10, 'LineWidth', 2);

% move the mobile device
while t < total_time
    % calculate the cell id of current position
    distances = sqrt((all_x - x).^2 + (all_y - y).^2);
    [~, cell_id] = min(distances); % ignore the first parameter(distance)
    
    % move the mobile device
    x_old = x;
    y_old = y;
    x = x + velocity * t_interval * cos(direction);
    y = y + velocity * t_interval * sin(direction);
    
    % check if the mobile device is still inside the map
    % =====================================================================
    if x < min(all_x) || x > max(all_x) || y < min(all_y) || y > max(all_y)
        % adjust the position and direction of the mobile device
        x = mod(x, 2*max(all_x)) - max(all_x);
        y = mod(y, 2*max(all_y)) - max(all_y);
        direction = direction + pi;
    end
    
    % calculate the cell id of new position
    distances = sqrt((all_x - x).^2 + (all_y - y).^2);
    [~, new_cell_id] = min(distances);
    
    % check if handoff occurs
    if new_cell_id ~= cell_id
        handoff_record = [handoff_record; t+ t_interval/2, cell_id, new_cell_id];
    end
    
    % plot the trajectory of mobile device
    plot([x_old, x], [y_old, y], 'o-', 'MarkerSize', 10, 'LineWidth', 2);
    
    % update time and parameters
    t = t + t_interval;
    velocity = 1 + 14*rand();
    direction = 2*pi*rand();
    t_interval = 1 + 5*rand();
end
hold off;

disp('  Time       Src_Cell  Dest_Cell');
disp(handoff_record);

% print the HO times
ho_count = size(handoff_record, 1);
disp(['Total handoff count: ', num2str(ho_count)]);

%================= functions to create the cluster =========================
%define the function of the basic 19 cells cluster
function plot_circles_and_labels(radius, all_x, all_y)
    grid on;
    for i = 1:19
        plot(all_x(i)+radius*cosd(0:60:360), all_y(i)+radius*sind(0:60:360), 'k');
        text(all_x(i), all_y(i)+100, num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        scatter(all_x(i), all_y(i), 'filled', 'MarkerFaceColor', 'b');
    end
    axis equal;
end

%----------------------for plot all cells-------------------------------
function drawCell(radius, x, y, label)
    plot(x+radius*cosd(0:60:360),y+radius*sind(0:60:360),'c');
    text(x, y+100, label, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    scatter(x, y, 'filled', 'MarkerFaceColor', 'b');
end

function plot_all(radius, isd, all_x, all_y)
    for i = 1:6
        offsetX = (sqrt((15*radius/2)^2 + 250^2)/isd)*isd*cosd(60*i-60-(atan(sqrt(3)/15)*180/pi));        offsetY = (sqrt((15*radius/2)^2 + 250^2)/isd)*isd*sind(60*i-60-(atan(sqrt(3)/15)*180/pi));
        for k = 1:19
            temp = 19*i+k;
            all_x(temp) = all_x(k)+offsetX;
            all_y(temp) = all_y(k)+offsetY;
            drawCell(radius, all_x(k)+offsetX, all_y(k)+offsetY, num2str(k));
        end
    end
    axis equal;
end