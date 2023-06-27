% Initialize arrays for distance and angle from center location to each cell site
numCellSites = 19;
siteDistances = zeros(1,numCellSites);
siteAngles = zeros(1,numCellSites);

isd = 500; % Inter-site distance
siteDistances(2:7) = isd;
siteAngles(2:7) = 30:60:360;

siteDistances(8:13) = 2*isd*cosd(30);
siteAngles(8:13) = 0:60:300;

siteDistances(14:19) = 2*isd;
siteAngles(14:19) = 30:60:360;

% Set center coordinates of the hexagonal cells
centerX = 0;
centerY = 0;

% Set radius of the cells
cellRadius = isd/sqrt(3);

% Initialize arrays for x and y coordinates of each cell site
xCoords = zeros(1,numCellSites*7);
yCoords = zeros(1,numCellSites*7);

% Calculate x and y coordinates
for i = 1:19
    xCoords(i) = centerX + siteDistances(i)*cosd(siteAngles(i));
    yCoords(i) = centerY + siteDistances(i)*sind(siteAngles(i));
end

% Create plot of the hexagonal cells
figure;
hold on;

% Plot original hexagonal pattern
for i = 1:19
    plot(xCoords(i)+cellRadius*cosd(0:60:360),yCoords(i)+cellRadius*sind(0:60:360),'k');
    text(xCoords(i)+100, yCoords(i), num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    scatter(xCoords(i), yCoords(i), 'filled', 'MarkerFaceColor', 'b');
end

% Loop for duplicating the figure and placing it around the current figure
for i = 1:6
    % Calculate the x and y offsets for each duplication
    offsetX = (sqrt((15*500/(2*sqrt(3)))^2+250^2)/500)*isd*cosd(60*i-60-(atan(sqrt(3)/15)*180/pi));
    offsetY = (sqrt((15*500/(2*sqrt(3)))^2+250^2)/500)*isd*sind(60*i-60-(atan(sqrt(3)/15)*180/pi));
    
    % Plot duplicated hexagonal pattern with offsets
    for k = 1:19
        temp = 19*i+k;
        xCoords(temp) = xCoords(k)+offsetX;
        yCoords(temp) = yCoords(k)+offsetY;
        plot(xCoords(k)+offsetX+cellRadius*cosd(0:60:360),yCoords(k)+offsetY+cellRadius*sind(0:60:360),'r');
        text(xCoords(k)+offsetX+100, yCoords(k)+offsetY, num2str(k), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        scatter(xCoords(k)+offsetX, yCoords(k)+offsetY, 'filled', 'MarkerFaceColor', 'b');
    end
end

axis equal;

%________________UPGRADED_CODE_ABOVE_______________%
%______________________SETUP_______________________%


BW = 10e6; % Channel bandwidth (Hz)
P_b = 33; % Power of base station (dBm)
P_m = 23; % Power of mobile device (dBm)
Gt = 14; % Transmitter antenna gain (dB)
Gr = 14; % Receiver antenna gain (dB)
ht_b = 50+1.5; % Height of base station (m)
ht_m = 1.5; % Height of mobile device (m)
T = 27+273.15; % Temperature (K)
k = 1.38e-23; % Boltzmann constant (J/K)

%Two Ray Ground Model --> g(d) = (ht*hr)^2 / d^4

% Define the parameters of the random walk mobility model
minDirection = 0;
maxDirection = 2*pi;
minSpeed = 1;
maxSpeed = 15; %velocity between [minSpeed, maxSpeed]
minT = 1; %mobile device moves t seconds
maxT = 6; 
totalTime = 900; %total sim time = 900 seconds
randSpeed = randi([minSpeed maxSpeed]);
randTime = randi([minT maxT]);
randDir = 2*pi*rand();
initialLocation = [250, 0]; % Initial location


%move direction from [0,2pi] 
%______________Variables_Above___________%
currentTime = 0;
currentLocation = initialLocation;
currentCell = checkCell(initialLocation(1),initialLocation(2));
handoffEvents = [];
handoff_amount = 0;
while currentTime < totalTime
    % Choose direction uniformly between [0, 2*pi]
    direction = rand() * 2 * pi;
    % Choose velocity uniformly between [minSpeed, maxSpeed]

    velocity = minSpeed + rand * (maxSpeed - minSpeed); 
    %velocity = randi([minSpeed maxSpeed])

    % Choose travel time uniformly between [minT, maxT]
    
    travelTime = minT + rand * (maxT - minT);
    %randTime = randi([minT maxT])

    % Update location based on direction, velocity, and travel time
    deltaX = velocity * cos(direction) * travelTime;
    deltaY = velocity * sin(direction) * travelTime;
    currentLocation = currentLocation + [deltaX, deltaY];
    
    newCell = checkCell(currentLocation(1), currentLocation(2));
    if newCell ~= currentCell
        handoff_amount = handoff_amount+1;
        % Record handoff event with current time, cell source, and cell destination
        handoffEvents = [handoffEvents; currentTime, currentCell, newCell];
        % Display current location
        disp([num2str(handoff_amount),') ','Time: ', num2str(currentTime), ', Location: (', num2str(currentLocation(1)), ', ', num2str(currentLocation(2)), ')', ', Source_Cell: ',num2str(currentCell),', Dest_Cell: ',num2str(newCell)]);
        % Update current cell
        currentCell = newCell;
        
    end

    % Update current time
    currentTime = currentTime + travelTime;
    
end

disp(['Number of Handoff: ', num2str(size(handoffEvents,1))]);


% Write a function that check the corresponding cell
function corresCellID = checkCell(x,y)
    numCellSites = 19*7;
    siteDistances = zeros(1,numCellSites);%19 column of 0s
    siteAngles = zeros(1,numCellSites); %19 column of 0s 
    isd = 500; % Inter-site distance
    siteDistances(2:7) = isd; %distance from 2-7 is 500
    siteAngles(2:7) = 30:60:360;
    siteDistances(8:13) = 2*isd*cosd(30);
    siteAngles(8:13) = 0:60:300;
    siteDistances(14:19) = 2*isd;
    siteAngles(14:19) = 30:60:360;
    centerX = 0;
    centerY = 0;
    cellRadius = isd/sqrt(3);
    xCoords = zeros(1,numCellSites);
    yCoords = zeros(1,numCellSites);
    for i = 1:19
        xCoords(i) = centerX + siteDistances(i)*cosd(siteAngles(i));
        yCoords(i) = centerY + siteDistances(i)*sind(siteAngles(i));
    end

    for i = 1:6
    % Calculate the x and y offsets for each duplication
        offsetX = (sqrt((15*500/(2*sqrt(3)))^2+250^2)/500)*isd*cosd(60*i-60-(atan(sqrt(3)/15)*180/pi));
        offsetY = (sqrt((15*500/(2*sqrt(3)))^2+250^2)/500)*isd*sind(60*i-60-(atan(sqrt(3)/15)*180/pi));
        for k = 1:19
            temp = 19*i+k;
            xCoords(temp) = xCoords(k)+offsetX;
            yCoords(temp) = yCoords(k)+offsetY;
        end
    end
    corresCellID = 100;
    for i = 1:19
    % Generate random points in each hexagon
        if inpolygon(x,y,xCoords(i)+cellRadius*cosd(0:60:360),yCoords(i)+cellRadius*sind(0:60:360))
            corresCellID = i;
        end
    end
    if corresCellID == 100
        %disp(['Out of 19-cell boundary, reaching outside one, Location: (', num2str(x), ', ' ,num2str(y), ')']);
        for i = 1:6
            for k = 1:19
                temp = 19*i+k;
                if inpolygon(x,y,xCoords(temp)+cellRadius*cosd(0:60:360),yCoords(temp)+cellRadius*sind(0:60:360))
                    corresCellID = k;
                end  
            end
        end
    end

end