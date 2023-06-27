clc; clear; close all;
% set parameters
n_bs = 19;          % number of base stations
isd = 500;          % inter-site distance_all (in meters)
L = isd/sqrt(3);    %length of the hexagon
freq = 2.4e9;       % carrier frequency (in Hz)
ptx = 33-30;        % transmit power of base stations (in dB)
pm = 23-30;         % transmit power of mobile devices (in dB)
g_bs = 14;          % antenna gain of base stations (in dB)
g_m = 14;           % antenna gain of mobile devices (in dB)
h_bs = 50;          % height of base stations (in meters)
h_ms = 1.5;         % height of mobile devices (in meters)
num_devices = 50;

% ------------------------------------1.1-----------------------------------------
% Define hexagonal central cell
hex_vertices = L*[cosd(60:60:360); sind(60:60:360)]';
hex_vertices = [hex_vertices; hex_vertices(1,:)]; % repeat first vertex to close shape

% center of hexagon
x_c = 0;
y_c = 0;

% hexagon's vertices
x = x_c + L*cosd(0:60:360);
y = y_c + L*sind(0:60:360);

% since I use the function and device_x and device_y need to be changed
global device_x; 
global device_y;

device_x = zeros(num_devices, 1);
device_y = zeros(num_devices, 1);

plot_points(L, num_devices, x, y, 0)
plot_cell(x, y)

% ------------------------------------1.2-----------------------------------------
num_cell = 19;
BW = 10e6;
p_bs_W = to_value(ptx);
p_m_W = to_value(pm);
gt_W = to_value(g_bs);
gr_W = to_value(g_m);

% cell_distance and the angles for the cells
cell_distance = zeros(num_cell);
angles = zeros(num_cell);

index1 = 2:7;
index2 = 8:13;
index3 = 14:19;
cell_distance(index1) = isd;
cell_distance(index2) = 2*isd*cosd(30);
cell_distance(index3) = 2*isd;
angles(index1) = 30:60:360;
angles(index2) = 0:60:300;
angles(index3) = 30:60:360;

% calculate the coordinates of cell's x and y
global cell_x;
global cell_y;

cell_x = zeros(num_cell, 1);
cell_y = zeros(num_cell, 1);

for i = 1:num_cell
    cell_x(i) = x_c + cell_distance(i)*cosd(angles(i));
    cell_y(i) = y_c + cell_distance(i)*sind(angles(i));
end

% calculate the distance_all from 50 devices to 19 cells
global distance_all;
distance_all = []; % would be 50x19 matrix

SINR = get_SINR(num_devices, num_cell, h_bs, h_ms, p_m_W, gt_W, gr_W); 
SC = get_shannon(BW, num_devices, SINR);

plot_shannon(distance_all, SC, 0)

% ------------------------------------1.3-----------------------------------------
time = 1000;   % total simulation time
buffer = 6e6;  % BS traffic buffer bits
CBR = [0.25e6, 0.5e6, 1e6]; % constant bit rate with [low, medium, high]

low_loss_rate = get_loss_rate(SC, CBR(1), buffer, time, 'constant', num_devices);
mid_loss_rate = get_loss_rate(SC,CBR(2), buffer, time, 'constant', num_devices);
high_loss_rate = get_loss_rate(SC,CBR(3), buffer, time, 'constant', num_devices);

plot_histogram(CBR, low_loss_rate, mid_loss_rate, high_loss_rate, 0)

% ------------------------------------B.1-----------------------------------------
plot_points(L, num_devices, x, y, 1)
plot_cell(x, y)

% ------------------------------------B.2-----------------------------------------
% Would change the original variable in 1-2 since I use the variable name
% calculate the Interference
SINR = get_SINR(num_devices, num_cell, h_bs, h_ms, p_m_W, gt_W, gr_W);
SC = get_shannon(BW, num_devices, SINR);

plot_shannon(distance_all, SC, 1)

% ------------------------------------B.3-----------------------------------------
lambda = CBR; % use the same bps as CBR

low_loss_rate = get_loss_rate(SC, lambda(1), buffer, time, 'poisson', num_devices);
mid_loss_rate = get_loss_rate(SC, lambda(2), buffer, time, 'poisson', num_devices);
high_loss_rate = get_loss_rate(SC,lambda(3), buffer, time, 'poisson', num_devices);

plot_histogram(CBR, low_loss_rate, mid_loss_rate, high_loss_rate, 1)

% ================= functions used in this HW =========================
function plot_points(L, num_devices, x, y, num)
    global device_x;
    global device_y;

    count = 0;
    while count < num_devices
        device_x_temp = rand * 2 * L - L;
        device_y_temp = rand * 2 * L - L;
        if inpolygon(device_x_temp, device_y_temp, x, y)
            count = count + 1;
            device_x(count) = device_x_temp;
            device_y(count) = device_y_temp;
        end
    end
    
    if num == 0
        figure('Name', 'Problem 1-1');
    else
        figure('Name', 'Problem B-1');
    end
    hold on;
    scatter(0, 0, 100, 'red', 's', 'filled');
    scatter(device_x, device_y, 50, 'blue', 'r', 'filled');
end

function plot_cell(x, y)
    plot([x x(1)], [y y(1)], 'k');
    xlabel('Distance (m)');
    ylabel('Distance (m)');
    legend('Central Base Station', 'Mobile Devices');
    title('Location of Central Base Station and Mobile Devices');
    axis square;
    grid on;
    hold off;
end

function plot_shannon(distance_all, SC, num)
    if num == 0
        figure('Name','Problem 1-2');
    else
        figure('Name','Problem B-2');
    end
    
    plot(distance_all(:,1), SC, 'o')
    xlabel('Distance (m)');
    ylabel('Shannon Capacity (bps)');
    title('Shannon Capacity vs Distance');
    grid on;
end

function plot_histogram(CBR, low_loss_rate, mid_loss_rate, high_loss_rate, num)
    x_values = [CBR(1), CBR(2), CBR(3)];                        % bit rate
    y_values = [low_loss_rate, mid_loss_rate, high_loss_rate];  % loss rate

    if num == 0
        figure('Name','Problem 1-3');
    else
        figure('Name','Problem B-3');
    end
    
    bar(x_values, y_values);
    
    for i = 1:length(x_values)
        text(x_values(i), y_values(i), num2str(y_values(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
    
    grid on
    xlabel('Bit rate (bps)');
    ylabel('Loss rate');
    title('Different Bit Rate vs Loss Rate');
    ylim([0, 1]);
end

function SINR = get_SINR(num_devices, num_cell, h_bs, h_ms, p_m_W, gt_W, gr_W)
    BW = 10e6;               % channel bandwidth (in Hz)
    T = 27+273.15;           % ambient temperature (in degree Kalvin)
    k = 1.38e-23;            % Boltzman's constant
    N = k*T*BW;              % noise value
    global device_x;
    global device_y;
    global cell_x;
    global cell_y;
    global distance_all;

    for i = 1:num_devices
       for j = 1:num_cell
            dx = device_x(i) - cell_x(j);
            dy = device_y(i) - cell_y(j);
            distance_all(i, j) = sqrt(dx^2 + dy^2);
        end
    end

    gd = ((h_bs*h_ms)^2)./distance_all.^4;
    Pr_W = gd.*p_m_W*gt_W*gr_W;

    I = zeros(size(Pr_W)); % interference
    for i = 1:size(Pr_W,1) 
       for j = 1:size(Pr_W,2)
            I(i,j) = sum(Pr_W(i, [1:j-1, j+1:end]));
        end
    end

    SINR = Pr_W./(I+N); % SINR in value
    return
end

function SC = get_shannon(BW, num_devices, SINR)
    each_BW = BW/num_devices;
    SC = zeros(50,1); % shannon capacity (bits/s); SC is a matrix(50x1)
   
    for i = 1:50
        SC(i,1) = each_BW*log2(1+SINR(i,1));
    end
    return
end

function loss_rate = get_loss_rate(SC, rate , buffer, time, type, num_devices)
    total_bit = 0;
    total_buffer_bit = 0;

    for i = 1:time
        % if the type is Poisson, we change the rate to poisson distribution
        if strcmp(type, 'poisson')
            rate = poissrnd(rate);
        end
        
        for i = 1:num_devices
            total_bit = total_bit + rate;
            max_capacity = SC(i,1);
            if rate > max_capacity 
                total_buffer_bit = total_buffer_bit + (rate-max_capacity);
            end
        end
    end

    loss_bit = total_buffer_bit - buffer;
    if loss_bit < 0
        loss_bit = 0;
    end

    loss_rate = loss_bit/total_bit;
    return
end


function result_value = to_value(db)
    result_value = 10^(db/10);
end