%------------------------1.1-------------------------
clc; clear; close all;
% set parameters
n_bs = 19;          % number of base stations
isd = 500;  % inter-site distance (in meters)

% hexagon's length
L = isd/sqrt(3);  %length of the hexagon

freq = 2.4e9;       % carrier frequency (in Hz)
bw = 10e6;          % channel bandwidth (in Hz)
ptx = 33-30;           % transmit power of base stations (in dB)
pm = 23-30;            % transmit power of mobile devices (in dB)

g_bs = 14;          % antenna gain of base stations (in dBi)
g_m = 14;           % antenna gain of mobile devices (in dBi)
h_bs = 50;         % height of base stations (in meters)
h_ms = 1.5;         % height of mobile devices (in meters)
T = 27+273.15; % ambient temperature (in degree Kalvin)
num_points = 50;
k = 1.38e-23;       % Boltzman's constant

%These three are in value
% ptx_w = (10^(ptx/10)/1000); 
% g_bs_w = 10^(g_bs/10);
% g_m_w = 10^(g_m/10);

% Define hexagonal central cell
hex_vertices = L*[cosd(60:60:360); sind(60:60:360)]';
hex_vertices = [hex_vertices; hex_vertices(1,:)]; % repeat first vertex to close shape

% For all the 19 cells
centers = [];
layer = 2;
center_locations = [0,0; 0, L*sqrt(3); L*sqrt(3)*sqrt(3)/2, L*sqrt(3)/2; L*sqrt(3)*sqrt(3)/2, -L*sqrt(3)/2;
                   0, -sqrt(3)*L; -L*sqrt(3)*sqrt(3)/2, -L*sqrt(3)/2; -L*sqrt(3)*sqrt(3)/2, L*sqrt(3)/2;
                   0,L*sqrt(3)*2; L*sqrt(3)*sqrt(3), L*sqrt(3); L*sqrt(3)*sqrt(3), -L*sqrt(3);
                   0, -sqrt(3)*L*2; -L*sqrt(3)*sqrt(3), -L*sqrt(3); -L*sqrt(3)*sqrt(3), L*sqrt(3);
                   L*sqrt(3)*sqrt(3)/2, L*sqrt(3)*3/2; L*3,0; L*sqrt(3)*sqrt(3)/2, -L*sqrt(3)*3/2;
                   -L*sqrt(3)*sqrt(3)/2, -L*sqrt(3)*3/2; -L*3,0; -L*sqrt(3)*sqrt(3)/2, L*sqrt(3)*3/2;]; 
%------------------------1.1-------------------------
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


%------------------------1.2-------------------------
% calculate distance between mobile devices and central BS
distance = sqrt((x_rand.^2) + (y_rand.^2));

% calculate path loss and trun into dB
gd = 10*log10(((h_bs*h_ms)^2) ./ (distance.^4));

% calculate received power at mobile devices
prx = ptx + g_bs + g_m + gd; 
% prx = ptx_w*g_bs_w*g_m_w.*gd;
% prx_dB = 10*log10(prx);

% plot the received power vs. distance
figure('Name', '1-2');
plot(distance, prx, 'o');
xlabel('Distance (m)');
ylabel('Received Power of MS(dB)');
title('Received Power of MS vs. Distance');
grid on;


%------------------------1.3-------------------------
I = 0; % caused by the other 18 cells(received power)
N = k * T * bw;
N_dB = 10*log10(N);

SINR = prx - (I+N_dB);

figure('Name','1-3');
scatter(distance, SINR, 'filled', 'MarkerFaceColor', 'b');
xlabel('Distance (m)');
ylabel('SINR (dB)');
title('SINR vs Distance');
grid on;


%------------------------2.1-------------------------
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
figure('Name', '2-1');
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


%------------------------2.2-------------------------
distance = sqrt((x_rand.^2) + (y_rand.^2));
gd = 10*log10(((h_bs*h_ms)^2) ./ (distance.^4));
prx_2 = pm + g_bs + g_m + gd; 

figure('Name', '2-2');
plot(distance, prx_2, 'o');
xlabel('Distance (m)');
ylabel('Received Power of BS (dB)');
title('Received Power of BS vs. Distance');
grid on;


%------------------------2.3-------------------------
SINR_2 = prx_2-(I+N_dB);

figure('Name','2-3');
scatter(distance, SINR_2, 'filled', 'MarkerFaceColor', 'b');
xlabel('Distance (m)');
ylabel('SINR (dB)');
title('SINR vs Distance');
grid on;


%------------------------B-1-------------------------
MS_position = [];
temp = num_points;
new_device_num = num_points;

while size(MS_position, 1) < num_points
    x = rand(temp, 1) *2*L-L;
    y = rand(temp, 1) *2*L-L;
    is_inside = inpolygon(x, y, hex_vertices(:,1), hex_vertices(:,2));
    MS_position = [MS_position; x(is_inside), y(is_inside)];
    temp = new_device_num-size(MS_position, 1);
end

for i = 2:19
    num_points = num_points+50;
    temp = 50;
    new_device_num = num_points;
    while size(MS_position, 1) < num_points
        x = rand(temp, 1) *2*L-L;
        y = rand(temp, 1) *2*L-L;
        is_inside = inpolygon(x, y, hex_vertices(:,1), hex_vertices(:,2));
        MS_position = [MS_position; x(is_inside)+center_locations(i,1), y(is_inside)+center_locations(i,2)];
        temp = new_device_num-size(MS_position, 1);
    end
end

for i = 1:size(center_locations,1)
    center = center_locations(i, :);
    start_idx = (i - 1) * 50 + 1;
    end_idx = i * 50;
    points = MS_position(start_idx:end_idx, :);
    diffs = points - center;
    distances(start_idx:end_idx) = sqrt(sum(diffs.^2, 2));
end

size(distances,1);
size(center_locations,1);

figure('Name','B-1');
scatter(center_locations(:,1), center_locations(:,2), 'filled', 'MarkerFaceColor', 'b');
hold on;
scatter(MS_position(:,1), MS_position(:,2), 'filled', 'MarkerFaceColor', 'r');
hold on;

plot(hex_vertices(:,1), hex_vertices(:,2), 'LineWidth', 2, 'Color', 'k');
hold on;
for i = 1:size(center_locations, 1)
    plot(hex_vertices(:,1)+center_locations(i,1), hex_vertices(:,2)+center_locations(i,2), 'LineWidth', 2, 'Color', 'k');
end

xlabel('x axis (m)');
ylabel('y axis (m)');
legend('Central BS', 'Mobile Devices', 'Other Cell');
title('Locations of 19 cells');
grid on;
hold off;


% ----------------------B-2-------------------------------
gd = 10*log10(((h_bs*h_ms)^2) ./ (distances.^4));
pr_dB_2 = pm + g_bs + g_m + gd; 

figure('Name','B-2');
scatter(distances,pr_dB_2, 'filled', 'MarkerFaceColor', 'r');
axis([0 max(distances+20) min(pr_dB_2-10) max(pr_dB_2+10)]);
xlabel('Distance (m)');
ylabel('Received Power (dB)');
title('Uplink Received Power v.s. Distance');
grid on
% ----------------------B-3-------------------------------
