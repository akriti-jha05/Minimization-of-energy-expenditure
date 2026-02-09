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

% Cray optimization algorithm Routing Function


function [Xbest, bestCost] = crayOptimizationAlgorithm(numSensorNodes, numGatewayNodes, maxIter, alpha, beta, baseStation, network_range)
    [network, nextHops, X, positions] = initializeNetwork(numSensorNodes, numGatewayNodes, network_range);
    numNodes = 60;
    dim = length(X); % Dimension of each position

    % Initialize crayfish
    crayfish(numNodes).position = [];
    crayfish(numNodes).bestPosition = [];
    crayfish(numNodes).bestCost = inf;
    Xshade = zeros(numNodes, dim); % Cave positions
    Xfood = rand(1, dim); % Food position

   globalBestPosition = rand(1, dim);
    globalBestCost = inf;

    % Initialize positions and best positions
    for i = 1:numNodes
        crayfish(i).position = initializeParticle(numSensorNodes, numGatewayNodes, nextHops, X, baseStation, positions);
        crayfish(i).bestPosition = rand(1, dim);
        crayfish(i).bestCost = fitnessFunction(crayfish(i).position, network, alpha, beta);

        % Update global best
        if crayfish(i).bestCost < globalBestCost
            globalBestCost = crayfish(i).bestCost;
            globalBestPosition = crayfish(i).bestPosition;
        end
    end

   

    % Main loop for COA optimization
    for iter = 1:maxIter
        % Temperature Definition
        C = 2 - (iter / maxIter); % Decreasing parameter
        temp = 20 + 15 * rand(); % Random temperature between 20 and 35
        
        for i = 1:numNodes
            Xnew = zeros(1, dim); % Initialize new position
            
            if temp > 30
                if rand < 0.5
                    % Summer Resort Stage
                    

                     Xshade(i, :) = (globalBestPosition + crayfish(i).bestPosition) / 2; % Update Xshade
                     for j = 1:dim 
                    % Ensure dimensional consistency
                       Xnew = crayfish(i).position(j) + C * rand() .* (Xshade(i, j) - crayfish(i).position(j)); % Eq.(6)
                     end
                else
                    % Competition Stage
                    z = randi(numNodes); % Randomly select another crayfish
                    Xshade(i, :) =  (globalBestPosition + crayfish(i).bestPosition); % Update Xshade
                    for j = 1:dim
                        Xnew(j) = crayfish(i).position(j) - crayfish(z).position(j) + Xshade(i, j); % Eq.(8)
                    end
                end
            else
                % Foraging Stage
                fitness_i = fitnessFunction(crayfish(i).position, network, alpha, beta);
                fitness_Xfood = fitnessFunction(Xfood, network, alpha, beta); % Fitness of the food position
                P = 3 * rand * fitness_i / fitness_Xfood; % Eq.(4) 
                
                if P > 2
                    % The food is too big
                    Xfood = exp(-1 / P) .* Xfood; % Eq.(12)
                    for j = 1:dim
                        Xnew(j) = crayfish(i).position(j) + cos(2 * pi * rand) * Xfood(j) * p_obj(temp) - sin(2 * pi * rand) * Xfood(j) * p_obj(temp); % Eq.(13)
                    end
                else
                   for j = 1:dim
                        Xnew(j) = (crayfish(i).position(j) - Xfood(j)) * p_obj(temp) + p_obj(temp) * rand() * crayfish(i).position(j); % Eq.(14)
                   end
                end
            end

            % Apply constraints to Xnew
            Xnew = max(0.00001, min(0.99999, Xnew));
            
            % Evaluate new position
            cost = fitnessFunction(Xnew, network, alpha, beta);

            % Update crayfish position
            if cost < crayfish(i).bestCost
                crayfish(i).bestCost = cost;
                crayfish(i).bestPosition = Xnew;
            end

            % Update global best
            if cost < globalBestCost
                globalBestCost = cost;
                globalBestPosition = Xnew;
            end

            % Update Xshade (cave position)
            Xshade(i, :) =  (globalBestPosition + crayfish(i).bestPosition);
        end
    end

    Xbest = globalBestPosition;
    bestCost = globalBestCost;
end

function y = p_obj(x)   %Eq.(4)
    y = 0.2*(1/(sqrt(2*pi)*3))*exp(-(x-25).^2/(2*3.^2));
end




%Fitness function to evaluate particle positions
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


%{
function cost = fitnessFunction(position, network, alpha, beta)
    maxDist = 0;
    maxHop = length(position) - 1;
    
    for i = 1:length(position) - 1
        if position(i) == position(i + 1)
            continue; % Skip consecutive same nodes
        end
        dist = network(position(i), position(i + 1));
        if dist > maxDist
            maxDist = dist;
        end
    end
    
    cost = alpha * maxDist + beta * maxHop;
end

%}

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
function [rounds, unusedGateways] = getRoundsUntilFirstNodeExhausted(numSensors, numGateways, baseStation, positions, E_elec, epsilon_fs, epsilon_mp, msg_size, d_0, initial_energy, nextHops, network, alpha, beta, Xbest)
    numNodes = numSensors + numGateways;
    energy_levels = ones(1, numNodes) * initial_energy;  % Initial energy levels for all nodes
    unusedGateways = zeros(1, numGateways);  % Track unused gateways
    
    rounds = 0;
    while all(energy_levels > 0)
        gatewaysUsedInRound = false(1, numGateways);  % Track which gateways are used in this round
        
        for sensor = 1:numSensors
            route = sensor;  % Start with the sensor node
            visited = false(1, numNodes);  % Track visited nodes to avoid loops
            visited(sensor) = true;
            
            while route(end) ~= baseStation
                current_node = route(end);
                validHops = nextHops{current_node};
                
                if isempty(validHops)
                    break;  % No valid hops available, terminate route
                end
                
                % Safeguard to ensure nextHopIdx does not exceed bounds
                X_idx = max(1, min(current_node, length(Xbest)));
                X_scaled = Xbest(X_idx);
                nextHopIdx = max(1, min(ceil(length(validHops) * X_scaled), length(validHops)));
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
                        gatewaysUsedInRound(node - numSensors) = true;  % Mark gateway as used
                        energy_levels(node) = energy_levels(node) - energy_cost;
                        
                        % Check if energy level drops to zero or below
                        if energy_levels(node) <= 0
                            unusedGateways = sum(~gatewaysUsedInRound);  % Count unused gateways in the final round
                            return;  % Stop if any node runs out of energy
                        end
                    end
                end
            end
        end
        
          % Count unused gateways for this round
        rounds = rounds + 1;
    end
    unusedGateways = sum(~gatewaysUsedInRound);
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

% Varying number of sensors from 100 to 700
sensor_counts = 100:100:700;
round_results_cray_optimizer = zeros(length(sensor_counts), size(base_station_locations, 1));

% Iterate over each base station location

for i = 1:size(base_station_locations, 1)
    for j = 1:length(sensor_counts)
        numSensors = sensor_counts(j);
        base_station = base_station_locations(i, :);
       
         [network, nextHops, X, positions] = initializeNetwork(numSensors, numGateways, network_range);
        % Initialize network and GWO routing
        [Xbest, bestCost] = crayOptimizationAlgorithm(numSensors, numGateways, numIterations, alpha, beta, base_station, network_range);
        
        % Calculate energy exhaustion rounds
        [rounds , unused_gateways] = getRoundsUntilFirstNodeExhausted(numSensors, numGateways, base_station, positions, E_elec, epsilon_fs, epsilon_mp, msg_size, d_0, initial_energy, nextHops, network, alpha, beta, Xbest);
        round_results_cray_optimizer(j,i) = rounds;
        % Store the result
        unused_gateways_results(j,i) = unused_gateways;
    end
end
% Plotting as a bar graph
%{
figure;
bar(sensor_counts, round_results);
xlabel('Number of Sensors');
ylabel('Number of Rounds before Energy Exhaustion');
title('Rounds vs. Sensors for Different Base Station Locations using COA');
legend('Base Station (500, 500)', 'Base Station (250, 0)', 'Base Station (400, 100)', 'Location', 'northwest');
grid on;

figure;
bar(sensor_counts, unused_gateways_results);
xlabel('Number of Sensors');
ylabel('Number of Unused Gateways');
title('Unused Gateways vs. Sensors for Different Base Station Locations using COA');
legend('Base Station (500, 500)', 'Base Station (250, 0)', 'Base Station (400, 100)', 'Location', 'northwest');
grid on;
%}
save('cray_optimizer.mat', 'round_results_cray_optimizer');