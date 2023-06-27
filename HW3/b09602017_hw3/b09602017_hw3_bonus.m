clc; clear; close all;
%------------------------1.1-------------------------
num_cell = 19;
isd = 500;
radius = isd/sqrt(3);

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

% calculate the coordinates of x and y
x_center = 250;
y_center = 0;

all_x = zeros(num_cell*7, 1);
all_y = zeros(num_cell*7, 1);

for i = 1:num_cell
    all_x(i) = x_center + distances(i)*cosd(angles(i));
    all_y(i) = y_center + distances(i)*sind(angles(i));
end

% create plot of the hexagonal cells
figure;
hold on;
plot_circles_and_labels(radius, all_x, all_y)
plot_all(radius, isd, all_x, all_y)

all_devices = 100;
all_points = [];
distance = [];

while size(all_points,1)<all_devices % do 100 times
    x = 1200 * (2*rand()-1);
    y = 1200 * (2*rand()-1);
    for i = 1:19
        if inpolygon(x,y,all_x(i)+radius*cosd(0:60:360),all_y(i)+radius*sind(0:60:360))
            all_points = [all_points; x y];
            text(x+20, y+20, num2str(size(all_points,1)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
    end
end 

title('Initial Devices Location')
scatter(all_points(:,1),all_points(:,2),'filled','r');
axis equal;

%------------------------1.2-------------------------
BW = 10e6; 
p_bs = 33; 
p_m = 23; 
gt = 14; 
gr = 14; 
h_b = 50+1.5;
h_m = 1.5; 
T = 27+273.15;
k = 1.38e-23; 
N = k*T*BW; 

p_bs_W = to_value(p_bs-30);
p_m_W = to_value(p_m-30);
gt_W = to_value(gt);
gr_W = to_value(gr);

% to remember all the distances in a 100x133 matrix, 
dx = repmat(all_points(:,1), 1, 133) - repmat(all_x', 100, 1);
dy = repmat(all_points(:,2), 1, 133) - repmat(all_y', 100, 1);
distance = sqrt(dx.^2 + dy.^2);

gd = ((h_b*h_m)^2)./distance.^2;
Pr_W = gd.*p_m_W*gt_W*gr_W;

[rows, cols] = size(Pr_W); % Determine the size of the pr_W matrix

% Create a matrix with same size as pr_W, but with zeros along the diagonal
Pr_W_no_diagonal = Pr_W;
for i = 1:min(rows, cols)
    Pr_W_no_diagonal(i,i) = 0;
end

% Compute the sum of each row except for the diagonal element
Interference = sum(Pr_W_no_diagonal, 2); 
SINR = Pr_W./(Interference+N);
SINR_dB = to_dB(SINR);

% parameters of the moving ones
min_v = 1;
max_v = 15;
min_t = 1;
max_t = 6; 
total_t = 900;
all_points;

current_t = 0;
current_location = all_points;
current_cell = [];
current_cell_pos = [];

current_cell = arrayfun(@(x) check(SINR_dB(x,:)), 1:100, 'UniformOutput', false);
current_cell = vertcat(current_cell{:});

ho_events = [];
ho_amount = 0;
new_cell = [];

while current_t < total_t
    direction = [];
    velocity = [];
    travel_t = [];
    delta_x = [];
    delta_y = [];
    
    travel_t = [travel_t; min_t + rand * (max_t - min_t)];

    for i = 1:100
        direction = [direction; rand() * 2 * pi];
        velocity = [velocity;min_v + rand * (max_v - min_v)]; 
        
        delta_x(i,1) = velocity(i,1) * cos(direction(i,1)) * travel_t;
        delta_y(i,1) = velocity(i,1) * sin(direction(i,1)) * travel_t;
    end 

    current_location = current_location + [delta_x delta_y];

    for i = 1:100
        for j = 1:133
            dx = current_location(i, 1) - all_x(j);
            dy = current_location(i, 2) - all_y(j);
            distance(i, j) = sqrt(dx^2 + dy^2);
        end
    end

    gd = ((h_b*h_m)^2) ./ distance.^2; 
    Pr_W = gd.*p_m_W*gt_W*gr_W; 

    mask = ones(size(Pr_W));
    mask(:,j) = 0;
    Interference = sum(Pr_W .* mask, 2);

    SINR = Pr_W ./ (Interference+N);
    SINR_dB = to_dB(SINR);

    for i = 1:100
        new_cell(i,:) = check(SINR_dB(i,:));

        if new_cell(i,:) ~= current_cell(i,:)
            ho_amount = ho_amount+1;
            ho_events = [ho_events; current_t, current_cell(i,:), new_cell(i,:)];
            disp([num2str(ho_amount),') ','Time: ', num2str(current_t), ', Source_Cell: ', num2str(current_cell(i,:)),', Dest_Cell: ',num2str(new_cell(i,:))]);
            current_cell(i,:) = new_cell(i,:); 
        end
    end
    current_t = current_t + travel_t; 
end

disp(['Number of Handoff: ', num2str(size(ho_events,1))]);

%======================= define the functions ============================
function plot_circles_and_labels(radius, all_x, all_y)
    grid on;
    for i = 1:19
        plot(all_x(i)+radius*cosd(0:60:360), all_y(i)+radius*sind(0:60:360), 'k');
        text(all_x(i), all_y(i)+75, num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
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
        offset_x = (sqrt((15*radius/2)^2 + 250^2)/isd)*isd*cosd(60*i-60-(atan(sqrt(3)/15)*180/pi));        
        offset_y = (sqrt((15*radius/2)^2 + 250^2)/isd)*isd*sind(60*i-60-(atan(sqrt(3)/15)*180/pi));
        for k = 1:19
            temp = 19*i+k;
            all_x(temp) = all_x(k)+offset_x;
            all_y(temp) = all_y(k)+offset_y;
            drawCell(radius, all_x(k)+offset_x, all_y(k)+offset_y, num2str(k));
        end
    end
    axis equal;
end
%-----------------------------------------------------------------------

function cell_id_out = check(vec)
    [~, maxIdx] = max(vec);
    mapping = mod(maxIdx-1, 19) + 1;
    mapping(mapping == 19) = 1;
    cell_id_out = mapping;
end

function result_dB = to_dB(value)
    result_dB = 10 * log10(value);
end

function result_value = to_value(db)
    result_value = 10^(db/10);
end