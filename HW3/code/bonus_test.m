clc; clear; close all;
%------------------------1.1-------------------------
num_cell = 19;
isd = 500;
radius = isd/sqrt(3);

% calculate the distances and the angles
distances = zeros(num_cell);
angles = zeros(num_cell);

index1 = 2:7;
index2 = 8:13;
index3 = 14:19;

distances(index1) = isd;
angles(index1) = 30:60:360;

distances(index2) = 2*isd*cosd(30);
angles(index2) = 0:60:300;

distances(index3) = 2*isd;
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
BW = 10e6; % Channel bandwidth (Hz)
P_b = 33; % Power of base station (dBm)
P_m = 23; % Power of mobile device (dBm)
Gt = 14; % Transmitter antenna gain (dB)
Gr = 14; % Receiver antenna gain (dB)
ht_b = 50+1.5; % Height of base station (m)
ht_m = 1.5; % Height of mobile device (m)
T = 27+273.15; % Temperature (K)
k = 1.38e-23; % Boltzmann constant (J/K)
N=k*T*BW; %noise
P_b_W = (10^((P_b-30)/10));
P_m_W = (10^((P_m-30)/10));
Gt_W = 10^(Gt/10);
Gr_W = 10^(Gr/10);

%distance is a 100x133 matrix, 
% for i = 1:100
%     for j = 1:133
%         dx = all_points(i, 1) - all_x(j);
%         dy = all_points(i, 2) - all_y(j);
%         distance(i, j) = sqrt(dx^2 + dy^2);
%     end
% end

dx = repmat(all_points(:,1), 1, 133) - repmat(all_x', 100, 1);
dy = repmat(all_points(:,2), 1, 133) - repmat(all_y', 100, 1);
distance = sqrt(dx.^2 + dy.^2);


gd = ((ht_b*ht_m)^2)./distance.^2; %also 100x133(100mobilephone, 133Basestation)
Pr_W = gd.*P_m_W*Gt_W*Gr_W; %also 100x133(100mobilephone, 133Basestation)

%______calculate interference --> add up all others except targeted_____%
% Interference = zeros(size(Pr_W)); % Initialize Interference to zeros
% 
% for i = 1:size(Pr_W,1) % Iterate over rows
%     for j = 1:size(Pr_W,2) % Iterate over columns
%         % Compute the sum of all elements in the ith row of Pr except for the jth column
%         Interference(i,j) = sum(Pr_W(i,[1:j-1,j+1:end]));
%     end
% end

% Determine the size of the Pr_W matrix
[rows, cols] = size(Pr_W);

% Create a matrix with same size as Pr_W, but with zeros along the diagonal
Pr_W_no_diagonal = Pr_W;
for i = 1:min(rows, cols)
    Pr_W_no_diagonal(i,i) = 0;
end

% Compute the sum of each row except for the diagonal element
Interference = sum(Pr_W_no_diagonal, 2);


%_________Interference_Calculated__________%
SINR = Pr_W./(Interference+N);
SINR_dB = 10*log10(SINR);

% Define the parameters of the random walk mobility model
minDirection = 0;
maxDirection = 2*pi;
minSpeed = 1;
maxSpeed = 15; %velocity between [minSpeed, maxSpeed]
minT = 1; %mobile device moves t seconds
maxT = 6; 
totalTime = 900; %total sim time = 900 seconds
% randSpeed = randi([minSpeed maxSpeed]);
% randTime = randi([minT maxT]);

randDir = 2*pi*rand();
randSpeed = 1 + rand()*14;
randTime = 1 + rand()*5;
all_points; % Initial locations

%______________Variables_Above___________%
currentTime = 0;
currentLocation = all_points;
currentCell = [];
currentCell_pos = [];

% for i = 1:100
%     currentCell = [currentCell; checkCell(SINR_dB(i,:))];
% end

currentCell = arrayfun(@(x) checkCell(SINR_dB(x,:)), 1:100, 'UniformOutput', false);
currentCell = vertcat(currentCell{:});


handoffEvents = [];
handoff_amount = 0;
newCell = [];
while currentTime < totalTime
    % Choose direction uniformly between [0, 2*pi]
    direction = [];
    velocity = [];
    travelTime = [];
    deltaX = [];
    deltaY = [];
    %parameters for position
    travelTime = [travelTime; minT + rand * (maxT - minT)];

    for i = 1:100
        direction = [direction; rand() * 2 * pi];
        velocity = [velocity;minSpeed + rand * (maxSpeed - minSpeed)]; 
        
        deltaX(i,1) = velocity(i,1) * cos(direction(i,1)) * travelTime;
        deltaY(i,1) = velocity(i,1) * sin(direction(i,1)) * travelTime;
    end 

    %update location
    currentLocation = currentLocation + [deltaX deltaY];

    %update distance from updated location
    for i = 1:100
        for j = 1:133
            dx = currentLocation(i, 1) - all_x(j);
            dy = currentLocation(i, 2) - all_y(j);
            distance(i, j) = sqrt(dx^2 + dy^2);
        end
    end

    %updated SINR
    gd = ((ht_b*ht_m)^2) ./ distance.^2; %also 100x133(100mobilephone, 133Basestation)
    Pr_W = gd.*P_m_W*Gt_W*Gr_W; %also 100x133(100mobilephone, 133Basestation)
    
    %______calculate interference --> add up all others except targeted_____%
%     Interference = zeros(size(Pr_W)); % Initialize Interference to zeros
%     for i = 1:size(Pr_W,1) % Iterate over rows
%         for j = 1:size(Pr_W,2) % Iterate over columns
%             % Compute the sum of all elements in the ith row of Pr except for the jth column
%             Interference(i,j) = sum(Pr_W(i,[1:j-1,j+1:end]));
%         end
%     end

    mask = ones(size(Pr_W));
    mask(:,j) = 0;
    Interference = sum(Pr_W .* mask, 2);


    %_________Interference_Calculated__________%
    SINR = Pr_W ./ (Interference+N);
    SINR_dB = 10 * log10(SINR);

    for i = 1:100
        newCell(i,:) = checkCell(SINR_dB(i,:));

        if newCell(i,:) ~= currentCell(i,:)
            handoff_amount = handoff_amount+1;
            handoffEvents = [handoffEvents; currentTime, currentCell(i,:), newCell(i,:)];
            disp([num2str(handoff_amount),') ','Time: ', num2str(currentTime), ', Source_Cell: ', num2str(currentCell(i,:)),', Dest_Cell: ',num2str(newCell(i,:))]);
            currentCell(i,:) = newCell(i,:); 
        end
    end
    currentTime = currentTime + travelTime; 
end

disp(['Number of Handoff: ', num2str(size(handoffEvents,1))]);

%define the function of the basic 19 cells cluster
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
        offsetX = (sqrt((15*radius/2)^2 + 250^2)/isd)*isd*cosd(60*i-60-(atan(sqrt(3)/15)*180/pi));
        offsetY = (sqrt((15*radius/2)^2 + 250^2)/isd)*isd*sind(60*i-60-(atan(sqrt(3)/15)*180/pi));
        for k = 1:19
            temp = 19*i+k;
            all_x(temp) = all_x(k)+offsetX;
            all_y(temp) = all_y(k)+offsetY;
            drawCell(radius, all_x(k)+offsetX, all_y(k)+offsetY, num2str(k));
        end
    end
    axis equal;
end

%----------------------for plot all cells-------------------------------

% function plot_all(radius, isd, all_x, all_y)
%     for i = 1:6
%     offsetX = (sqrt((15*500/(2*sqrt(3)))^2+250^2)/500)*isd*cosd(60*i-60-(atan(sqrt(3)/15)*180/pi));
%     offsetY = (sqrt((15*500/(2*sqrt(3)))^2+250^2)/500)*isd*sind(60*i-60-(atan(sqrt(3)/15)*180/pi));
%     for k = 1:19
%         temp = 19*i+k;
%         all_x(temp) = all_x(k)+offsetX;
%         all_y(temp) = all_y(k)+offsetY;
%         plot(all_x(k)+offsetX+radius*cosd(0:60:360),all_y(k)+offsetY+radius*sind(0:60:360),'c');
%         text(all_x(k)+offsetX+100, all_y(k)+offsetY, num2str(k), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%         scatter(all_x(k)+offsetX, all_y(k)+offsetY, 'filled', 'MarkerFaceColor', 'b');
%     end
%     end
% end

function corresCellID = checkCell(vec)
    [~, maxIdx] = max(vec);
    mapping = mod(maxIdx-1, 19) + 1;
    mapping(mapping == 19) = 1;
    corresCellID = mapping;
end

% function corresCellID = checkCell(vec)
%     [~, corresCellID] = max(vec);
%     temporary = corresCellID;
%     mapping = mod(temporary-1, 19) + 1;
%     if mapping == 19
%         mapping = 1;
%     end
%     corresCellID = mapping;
% end