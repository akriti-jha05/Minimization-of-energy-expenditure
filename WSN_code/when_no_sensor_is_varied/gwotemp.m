clc;
clear all;
close all;


network_range = 500;  % Range of the network
bandwidth = 1e6;  % 1 MB/s
msg_size = 525;   % 500 bytes + 25-byte header
E_elec = 50e-9;  % Energy per bit for electronics
epsilon_fs = 10e-12;  % Free space energy loss coefficient
epsilon_mp = 0.0013e-12;  % Multipath energy loss coefficient
d_0 = 75;
initial_energy = 10;

% Function to initialize the network
function [network, nextHops, X, positions] = initializeNetwork(numSensorNodes, numGatewayNodes, network_range)
    numNodes = numSensorNodes + numGatewayNodes;
    positions = rand(numNodes, 2) * network_range; % Random positions within the given range

    % Network adjacency matrix
    network = zeros(numNodes, numNodes);

    % Define range limit for sensor-to-gateway and gateway-to-gateway
    range_limit = 100;

    % Fill the network adjacency matrix
    for i = 1:numSensorNodes
        for j = numSensorNodes + 1:numNodes
            d = sqrt((positions(i, 1) - positions(j, 1))^2 + (positions(i, 2) - positions(j, 2))^2);
            if d < range_limit
                network(i, j) = d; % Sensor to Gateway
            end
        end
    end

    for i = numSensorNodes + 1:numNodes
        for j = numSensorNodes + 1:numNodes
            if i ~= j
                d = sqrt((positions(i, 1) - positions(j, 1))^2 + (positions(i, 2) - positions(j, 2))^2);
                if d < range_limit
                    network(i, j) = d; % Gateway to Gateway
                end
            end
        end
    end

    % Next hops for each node
    nextHops = cell(numNodes, 1);
    for i = 1:numNodes
        nextHops{i} = find(network(i, :) > 0);
    end

    X = rand(1, numNodes); % Random selection probabilities
end

% Function to initialize particle positions
function position = initializeParticle(numSensorNodes, numGatewayNodes, nextHops, X, baseStation, positions)
    numNodes = numSensorNodes + numGatewayNodes;
    position = zeros(1, numNodes + 1); % Initialize position vector

    for i = 1:numNodes
        if isempty(nextHops{i})
            position(i) = i; % Set node to itself if no next hops
        else
            nextHopIdx = ceil(length(nextHops{i}) * X(i));
            position(i) = nextHops{i}(min(nextHopIdx, length(nextHops{i}))); % Choose next hop
        end
    end

    % Find the nearest gateway to the base station coordinates
    gateway_positions = positions(numSensorNodes + 1:end, :);
    distances = sqrt((gateway_positions(:, 1) - baseStation(1)).^2 + (gateway_positions(:, 2) - baseStation(2)).^2);
    [~, baseStationIndex] = min(distances);
    baseStationIndex = baseStationIndex + numSensorNodes; % Adjust index for gateway nodes

    % Set the last element of position to the index of the nearest gateway to the base station
    position(end) = baseStationIndex;
end

% Grey Wolf Optimization (GWO) Routing Function
function [Xbest, bestCost] = gwoRouting(numSensorNodes, numGatewayNodes, maxIter, alpha, beta, baseStation, network_range)
    [network, nextHops, X, positions] = initializeNetwork(numSensorNodes, numGatewayNodes, network_range);
    numNodes = numSensorNodes + numGatewayNodes;

    % Initialize wolves (particles in PSO terms)
    wolves(numNodes).position = [];
    wolves(numNodes).bestPosition = [];
    wolves(numNodes).bestCost = inf;

    globalBestPosition = [];
    globalBestCost = inf;

    % Initialize positions and best positions
    for i = 1:numNodes
        wolves(i).position = initializeParticle(numSensorNodes, numGatewayNodes, nextHops, X, baseStation, positions);
        wolves(i).bestPosition = wolves(i).position;
        wolves(i).bestCost = fitnessFunction(wolves(i).position, network, alpha, beta);

        % Update global best
        if wolves(i).bestCost < globalBestCost
            globalBestCost = wolves(i).bestCost;
            globalBestPosition = wolves(i).bestPosition;
        end
    end

    % GWO parameters
    a = 2;  % Alpha parameter
    a_damp = 0.98;  % Alpha damping coefficient

    % Main loop for GWO optimization
    iter = 1;
    while iter <= maxIter
        % Determine alpha, beta, and delta positions
        [~, sortedIndices] = sort([wolves.bestCost]);
        alphaPosition = wolves(sortedIndices(1)).position;
        betaPosition = wolves(sortedIndices(2)).position;
        deltaPosition = wolves(sortedIndices(3)).position;

        for i = 1:numNodes
            % Update positions using GWO formula
            r1 = rand(1);
            r2 = rand();
            A1 = 2 * a * r1 - a;
            C1 = 2 * r2;
            D_alpha = abs(C1 * alphaPosition - wolves(i).position);
            X1 = alphaPosition - A1 * D_alpha;

            r1 = rand();
            r2 = rand();
            A2 = 2 * a * r1 - a;
            C2 = 2 * r2;
            D_beta = abs(C2 * betaPosition - wolves(i).position);
            X2 = betaPosition - A2 * D_beta;

            r1 = rand();
            r2 = rand();
            A3 = 2 * a * r1 - a;
            C3 = 2 * r2;
            D_delta = abs(C3 * deltaPosition - wolves(i).position);
            X3 = deltaPosition - A3 * D_delta;

            newPosition = (X1 + X2 + X3) / 3;

            % Apply constraints
            %wolves(i).position = max(0, min(1, wolves(i).position));
           if newPosition <= 0
        % Replace with a small random number tending to zero
                wolves(i).position = rand() * 0.001;  % Adjust 0.01 as needed
           elseif newPosition >= 1
                wolves(i).position = 1;
           else
                wolves(i).position = newPosition;
           end
            % Evaluate cost function
            cost = fitnessFunction(wolves(i).position, network, alpha, beta);

            % Update best position and global best
            if cost < wolves(i).bestCost
                wolves(i).bestCost = cost;
                wolves(i).bestPosition = wolves(i).position;
            end

            if cost < globalBestCost
                globalBestCost = cost;
                globalBestPosition = wolves(i).bestPosition;
            end
        end

        % Update alpha parameter
        a = 2 -iter*(2/maxIter);
        iter = iter + 1;
    end

    Xbest = globalBestPosition;
    bestCost = globalBestCost;
end

% Fitness function to evaluate particle positions
function cost = fitnessFunction(position, positions, alpha, beta)
    maxDist = 0;
    maxHop = length(position) - 1;
    
    % Initialize route array (if needed)
    % route = zeros(1, length(position));
    
    numPositions = size(position, 1); % Get the number of positions
    
    for i = 1:(length(position) - 1)
        % Check if indices are within bounds
        if position(i) < 1 || position(i) > numPositions || position(i + 1) < 1 || position(i + 1) > numPositions
            % Skip invalid indices
            continue;
        end
        
        % Check if consecutive nodes are the same
        if position(i) == position(i + 1)
            continue; % Skip consecutive same nodes
        end
        
        % Calculate Euclidean distance between nodes
        node1 = position(i);
        node2 = position(i + 1);
        dist = sqrt((position(node2, 1) - position(node1, 1))^2 + (position(node2, 2) - position(node1, 2))^2);
        
        % Update maximum distance
        if dist > maxDist
            maxDist = dist;
        end
    end
    
    % Calculate cost based on max distance and hop count
    cost = alpha * maxDist + beta * maxHop;
end




% Function to calculate the energy cost of a route
function energy_cost = calculateEnergyCost(route, E_elec, epsilon_fs, epsilon_mp, msg_size, d_0, positions)
    numNodes = length(route);
    energy_cost = 0;

    for i = 1:numNodes - 1
        node1 = route(i);
        node2 = route(i + 1);
        
        if node1 == node2
            continue; % Skip if the route has consecutive same nodes
        end
        
        % Check if nodes are within bounds
        if node1 > size(positions, 1) || node2 > size(positions, 1)
            error('Invalid node index in route');
        end
        
        d = sqrt((positions(node2, 1) - positions(node1, 1))^2 + (positions(node2, 2) - positions(node1, 2))^2);
        
        if d < d_0
            energy_transmitter = E_elec * msg_size + epsilon_fs * (d^2) * msg_size;
        else
            energy_transmitter = E_elec * msg_size + epsilon_mp * (d^4) * msg_size;
        end
        
        energy_receiver = E_elec * msg_size;
        energy_cost = energy_cost + energy_transmitter + energy_receiver;
    end
end

% Function to calculate the number of rounds until the first node is exhausted
function rounds = getRoundsUntilFirstNodeExhausted(numSensors, numGateways, baseStation, positions, E_elec, epsilon_fs, epsilon_mp, msg_size, d_0, initial_energy, nextHops, network, alpha, beta, Xbest)
    numNodes = numSensors + numGateways;
    energy_levels = ones(1, numNodes) * initial_energy;  % Initial energy levels for all nodes
    
    rounds = 0;
    while all(energy_levels > 0)
        for sensor = 1:numSensors
            route = sensor;  % Start with the sensor node
            visited = false(1, numNodes);  % Track visited nodes to avoid loops
            visited(sensor) = true;
            
            while route(end) ~= baseStation
                current_node = route(end);
                if isempty(nextHops{current_node})
                    break;  % No next hops available, terminate route
                end
                
                % Safeguard to ensure nextHopIdx does not exceed bounds
                validHops = nextHops{current_node};
                if isempty(validHops)
                    break; % No valid hops, terminate route
                end
                
                % Determine the index for Xbest
                if current_node <= length(Xbest) && current_node >= 1
                    X_idx = current_node;  % Use current_node as index for Xbest
                else
                    X_idx = ceil(rand * length(Xbest));  % Random index if current_node is out of bounds
                end
                
                % Compute nextHopIdx using Xbest scaled to validHops length
                X_scaled = Xbest(X_idx);
                nextHopIdx = max(1, ceil(length(validHops) * X_scaled));
                nextHopIdx = min(nextHopIdx, length(validHops));  % Ensure nextHopIdx is within bounds
                next_node = validHops(nextHopIdx);

                if visited(next_node)
                    break;  % Avoid loops by breaking if the next node has already been visited
                end
                
                route = [route, next_node];  % Add the next node to the route
                visited(next_node) = true;
                
                % If the route is excessively long, break to avoid infinite loops
                if length(route) > numNodes
                    break;
                end
            end
            
            % Calculate energy cost for the current route
            if length(route) > 1
                energy_cost = calculateEnergyCost(route, E_elec, epsilon_fs, epsilon_mp, msg_size, d_0, positions);
                
                % Update energy levels for gateways in the route
                for node = route
                    if node > numSensors && node <= numNodes
                        % Deduct energy cost for gateways only
                        energy_levels(node) = energy_levels(node) - energy_cost;
                        
                        % Check if energy level drops to zero or below
                        if energy_levels(node) <= 0
                            return;  % Stop if any node runs out of energy
                        end
                    elseif node > numNodes
                        error('Invalid node index in route');
                    end
                end
            end
        end
        rounds = rounds + 1;
    end
end
searchagent = 30;
% Simulation Parameters
numIterations = 700;        % Maximum number of iterations
alpha = 0.5;          % Weight for distance
beta = 1 - alpha;     % Weight for hop count

base_station_locations = [
    500, 500;  % Scenario 1
      250,250
];

numGateways = 60;  % Number of gateways

% Varying number of sensors from 140 to 700
sensor_counts = 100:100:700;
round_results_gwo = zeros(length(sensor_counts), size(base_station_locations, 1));

% Iterate over each base station location
for i = 1:size(base_station_locations, 1)
    for j = 1:length(sensor_counts)
        numSensors = sensor_counts(j);
        
        base_station = base_station_locations(i, :);
         [network, nextHops, X, positions] = initializeNetwork(numSensors, numGateways, network_range);
        % Initialize network and GWO routing
        [Xbest, bestCost] = gwoRouting(numSensors, numGateways, numIterations, alpha, beta, base_station, network_range);
        
        % Calculate energy exhaustion rounds
        rounds = getRoundsUntilFirstNodeExhausted(numSensors, numGateways, base_station, positions, E_elec, epsilon_fs, epsilon_mp, msg_size, d_0, initial_energy, nextHops, network, alpha, beta, Xbest);
        round_results_gwo(j, i) = rounds;
    end
end

% Plotting as a bar graph
%{
figure;
bar(sensor_counts, round_results);
xlabel('Number of Sensors');
ylabel('Number of Rounds before Energy Exhaustion');
title('Rounds vs. Sensors for Different Base Station Locations using GWO');
legend('Base Station (500, 500)', 'Base Station (250, 0)', 'Base Station (400, 100)', 'Location', 'northwest');
grid on;
%}
save('grey_wolf.mat', 'round_results_gwo');