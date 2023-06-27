ISD = 500; % Inter site distance (m)
BW = 10e6; % Channel bandwidth (Hz)
P_b = 33; % Power of base station (dBm)
P_m = 23; % Power of mobile device (dBm)
Gt = 14; % Transmitter antenna gain (dB)
Gr = 14; % Receiver antenna gain (dB)
ht_b = 50+1.5; % Height of base station (m)
ht_m = 1.5; % Height of mobile device (m)
T = 27+273.15; % Temperature (K)
k = 1.38e-23; % Boltzmann constant (J/K)

% Calculate the thermal noise power
N = k*T*BW; % Thermal Noise Power  (W/Hz)


%Q1-1
mobile_device_num = 50; % number of mobile devices
cell_radius = (250*2)/sqrt(3); % inter site distance
central_bs_location = [0, 0]; % coordinates of central BS

% Define hexagonal central cell
hex_vertices = cell_radius*[cosd(60:60:360); sind(60:60:360)]';
hex_vertices = [hex_vertices; hex_vertices(1,:)]; % repeat first vertex to close shape

% Generate random positions for mobile devices in central cell
mobile_device_positions = [];
temp = mobile_device_num;
while size(mobile_device_positions, 1) < mobile_device_num
    % Generate random positions in bounding box of hexagonal cell
    x = (rand(temp, 1) - 0.5) * cell_radius*2;
    y = (rand(temp, 1) - 0.5) * cell_radius*2;
    % Check if point is inside hexagonal cell
    is_inside = inpolygon(x, y, hex_vertices(:,1), hex_vertices(:,2));
    mobile_device_positions = [mobile_device_positions; x(is_inside), y(is_inside)];
    temp = 50-size(mobile_device_positions, 1);
end
distances = sqrt(sum((mobile_device_positions - central_bs_location).^2, 2));

% Plot locations of central BS, hexagonal cell, and mobile devices
figure('Name','Question 1-1');
scatter(central_bs_location(1), central_bs_location(2), 'filled', 'MarkerFaceColor', 'b');
hold on;
scatter(mobile_device_positions(:,1), mobile_device_positions(:,2), 'filled', 'MarkerFaceColor', 'r');
plot(hex_vertices(:,1), hex_vertices(:,2), 'LineWidth', 2, 'Color', 'k');
axis([-cell_radius-15 cell_radius+15 -cell_radius-15 cell_radius+15]);
xlabel('x-axis (m)');
ylabel('y-axis (m)');
legend('Central BS', 'Mobile Devices', 'Hexagonal Cell');
title('Locations of Central BS, Hexagonal Cell, and Mobile Devices');

%Q1-2 y-axis: received power x-axis:distance

%received power = g(d)*P_b*Gt*Gr where g(d) = ((ht_b)*(ht_m))^2/d^4
gd = ((ht_b*ht_m)^2)./(distances.^4);

%Pr_dB = P_b+Gt+Gr-gd
P_b_W = (10^(P_b/10)/1000);
Gt_W = 10^(Gt/10);
Gr_W = 10^(Gr/10);
Pr = P_b_W*Gt_W*Gr_W.*gd;
Pr_dB = 10*log10(Pr);

figure('Name','Question 1-2');
scatter(distances, Pr_dB, 'filled', 'MarkerFaceColor', 'r');
axis([0 max(distances+100) min(Pr_dB-10) max(Pr_dB+10)]);
xlabel('Distance (m)');
ylabel('Received Power (dB)');
title('Downlink Received Power vs Distance');

%Q1-3
%SINR = S/(I+N) N - Themal noise, S - Received Power 
I = 0;
SINR = Pr./(I+N);
SINR_dB = 10*log10(SINR);
figure('Name','Question 1-3');
scatter(distances, SINR_dB, 'filled', 'MarkerFaceColor', 'b');
axis([0 max(distances+100) min(SINR_dB-10) max(SINR_dB+10)]);
xlabel('Distances (m)');
ylabel('SINR (dB)');
title('SINR vs Distance');


%Q2-1
figure('Name','Question 2-1');
scatter(central_bs_location(1), central_bs_location(2), 'filled', 'MarkerFaceColor', 'b');
hold on;
scatter(mobile_device_positions(:,1), mobile_device_positions(:,2), 'filled', 'MarkerFaceColor', 'r');
plot(hex_vertices(:,1), hex_vertices(:,2), 'LineWidth', 2, 'Color', 'k');
axis([-cell_radius-15 cell_radius+15 -cell_radius-15 cell_radius+15]);
xlabel('x-axis (m)');
ylabel('y-axis (m)');
legend('Central BS', 'Mobile Devices', 'Hexagonal Cell');
title('Locations of Central BS and Mobile Devices');


%Q2-2
P_m_W = (10^(P_m/10)/1000);
Gt_W = 10^(Gt/10);
Gr_W = 10^(Gr/10);
Pr_2 = P_m_W*Gt_W*Gr_W.*gd;
Pr_dB_2 = 10*log10(Pr_2);

figure('Name','Question 2-2');
scatter(distances, Pr_dB_2, 'filled', 'MarkerFaceColor', 'r');
axis([0 max(distances+100) min(Pr_dB_2-10) max(Pr_dB_2+10)]);
xlabel('Distance (m)');
ylabel('Received Power (dB)');
title('Uplink Received Power vs Distance');


gd = ((ht_b*ht_m)^2)./(distances.^4);
P_b_W = (10^(P_b/10)/1000);
Gt_W = 10^(Gt/10);
Gr_W = 10^(Gr/10);
Pr = P_b_W*Gt_W*Gr_W.*gd;
Pr_dB = 10*log10(Pr);

figure('Name','Question 1-2');
scatter(distances, Pr_dB, 'filled', 'MarkerFaceColor', 'r');
axis([0 max(distances+100) min(Pr_dB-10) max(Pr_dB+10)]);
xlabel('Distance (m)');
ylabel('Received Power (dB)');
title('Downlink Received Power vs Distance');


%Q2-3
SINR_2 = Pr_2./(I+N);
SINR_dB_2 = 10*log10(SINR_2);
figure('Name','Question 2-3');
scatter(distances, SINR_dB_2, 'filled', 'MarkerFaceColor', 'b');
axis([0 max(distances+100) min(SINR_dB_2-10) max(SINR_dB_2+10)]);
xlabel('Distances (m)');
ylabel('SINR (dB)');
title('Uplink SINR vs Distance');

%Bonus B-1