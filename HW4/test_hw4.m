clc; clear; close all;
% set parameters
n_bs = 19;          % number of base stations
isd = 500;          % inter-site distance (in meters)
L = isd/sqrt(3);    %length of the hexagon
freq = 2.4e9;       % carrier frequency (in Hz)
ptx = 33-30;        % transmit power of base stations (in dB)
pm = 23-30;         % transmit power of mobile devices (in dB)
g_bs = 14;          % antenna gain of base stations (in dB)
g_m = 14;           % antenna gain of mobile devices (in dB)
h_bs = 50;          % height of base stations (in meters)
h_ms = 1.5;         % height of mobile devices (in meters)
num_points = 50;

% ------------------------------------1.1-----------------------------------------
% Define hexagonal central cell
hex_vertices = L*[cosd(60:60:360); sind(60:60:360)]';
hex_vertices = [hex_vertices; hex_vertices(1,:)]; % repeat first vertex to close shape

% center of hexagon
x_c = 0;
y_c = 0;

% center cell's vertices
x = x_c + L*cosd(0:60:360);
y = y_c + L*sind(0:60:360);

% generate 50 points that's in the cell 50x1
device_x = zeros(num_points, 1);
device_y = zeros(num_points, 1);

% calculate the 50 devices' coordinations
count = 0;
while count < num_points
    device_x_temp = rand * 2 * L - L;
    device_y_temp = rand * 2 * L - L;
    if inpolygon(device_x_temp, device_y_temp, x, y)
        count = count + 1;
        device_x(count) = device_x_temp;
        device_y(count) = device_y_temp;
    end
end

% plot the points in the hexagon
figure('Name', '1-1');
hold on;
scatter(0, 0, 100, 'red', 's', 'filled');  % plot the central base station
scatter(device_x, device_y, 50, 'blue', 'r', 'filled');

% plot hexagon
plot([x x(1)], [y y(1)], 'k');
xlabel('Distance (m)');
ylabel('Distance (m)');
legend('Central Base Station', 'Mobile Devices', 'Cell Region');
title('Location of Central Base Station and Mobile Devices');
axis square;
grid on;
hold off;

% ------------------------------------1.2-----------------------------------------
% parameters for SINR
BW = 10e6;               % channel bandwidth (in Hz)
each_BW = BW/50;
T = 27+273.15;           %  ambient temperature (in degree Kalvin)
k = 1.38e-23;            % Boltzman's constant
N = k*T*BW; % noise value

num_cell = 19;
p_bs_W = to_value(ptx);
p_m_W = to_value(pm);
gt_W = to_value(g_bs);
gr_W = to_value(g_m);

% remember the distances and the angles
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

% calculate the coordinates of cell's x and y
cell_x = zeros(num_cell, 1);
cell_y = zeros(num_cell, 1);

for i = 1:num_cell
    cell_x(i) = x_c + distances(i)*cosd(angles(i));
    cell_y(i) = y_c + distances(i)*sind(angles(i));
end

% calculate the distance from 50 devices to 19 cells
distance = [];

for i = 1:50
    for j = 1:19
        dx = device_x(i) - cell_x(j);
        dy = device_y(i) - cell_y(j);
        distance(i, j) = sqrt(dx^2 + dy^2);
    end
end

gd = ((h_bs*h_ms)^2)./distance.^4;
Pr_W = gd.*p_m_W*gt_W*gr_W;

% calculate the interference
Interference = zeros(size(Pr_W));
for i = 1:size(Pr_W,1) 
    for j = 1:size(Pr_W,2)
        % Compute the sum of all elements in the ith row of Pr except for the jth column
        Interference(i,j) = sum(Pr_W(i,[1:j-1,j+1:end]));
    end
end

SINR = Pr_W./(Interference+N); % SINR in value

SC = zeros(50,1); % shannon capacity (bits/s)
for i = 1:50
    SC(i,1) = each_BW*log2(1+SINR(i,1));
end

figure('Name','Problem 1-2');
% scatter(distance(:,1), SC, 'filled', 'MarkerFaceColor', 'b');
plot(distance(:,1), SC, 'o')
xlabel('Distance (m)');
ylabel('Shannon Capacity (bps)');
title('Shannon Capacity vs Distance');
grid on;


function result_dB = to_dB(value)
    result_dB = 10*log10(value);
end

function result_value = to_value(db)
    result_value = 10^(db/10);
end

% ================= functions used in this HW =========================
% function plot_points(L, num_points, x, y, device_x, device_y)
%     count = 0;
%     while count < num_points
%         device_x_temp = rand * 2 * L - L;
%         device_y_temp = rand * 2 * L - L;
%         if inpolygon(device_x_temp, device_y_temp, x, y)
%             count = count + 1;
%             device_x(count) = device_x_temp;
%             device_y(count) = device_y_temp;
%         end
%     
%     end
%     
%     % plot the points in the hexagon
%     figure('Name', 'Problem 1-1');
%     hold on;
%     scatter(0, 0, 100, 'red', 's', 'filled');  % plot the central base station
%     scatter(device_x, device_y, 50, 'blue', 'r', 'filled');
% end
% 
% function plot_hexa(x, y)
%     plot([x x(1)], [y y(1)], 'k');
%     xlabel('Distance (m)');
%     ylabel('Distance (m)');
%     legend('Central Base Station', 'Mobile Devices');
%     title('Location of Central Base Station and Mobile Devices');
%     axis square;
%     grid on;
%     hold off;
% end