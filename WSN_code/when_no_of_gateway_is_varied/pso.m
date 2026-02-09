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

function position = initializeParticle(numSensorNodes, numGatewayNodes, nextHops, X, baseStation, positions)
    numNodes = numSensorNodes + numGatewayNodes;
    position = zeros(1, numSensorNodes + 1); % Initialize position vector

    for i = 1:numSensorNodes
        if isempty(nextHops{i})
            position(i) = i; % Set node to itself if no next hops
        else
            nextHopIdx = ceil(length(nextHops{i}) * X(i));
            position(i) = nextHops{i}(nextHopIdx); %  choose next hop
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

% Main PSO Routing Function
function [Xbest, bestCost] = psoRouting(numSensorNodes, numGatewayNodes, numParticles, maxIter, alpha, beta, baseStation, network_range)
    [network, nextHops, X, positions] = initializeNetwork(numSensorNodes, numGatewayNodes, network_range);
    numNodes = numSensorNodes + numGatewayNodes;

    % Initialize particles
    particles(numParticles).position = [];
    particles(numParticles).velocity = [];
    particles(numParticles).bestPosition = [];
    particles(numParticles).bestCost = inf;

    globalBestPosition = [];
    globalBestCost = inf;

    for i = 1:numParticles
        particles(i).position = initializeParticle(numSensorNodes, numGatewayNodes, nextHops, X, baseStation, positions);

        particles(i).velocity = zeros(1, numSensorNodes + 1); % Include the base station
        particles(i).bestPosition = particles(i).position;
        particles(i).bestCost = fitnessFunction(particles(i).position, network, alpha, beta);

        if particles(i).bestCost < globalBestCost
            globalBestCost = particles(i).bestCost;
            globalBestPosition = particles(i).bestPosition;
        end
    end

    w = 0.5; % Inertia weight
    c1 = 1.2;  % Cognitive (personal) weight
    c2 = 2;  % Social (global) weight

    for iter = 1:maxIter
        for i = 1:numParticles
            inertia = w * particles(i).velocity;
            cognitive = c1 * rand * (particles(i).bestPosition - particles(i).position);
            social = c2 * rand * (globalBestPosition - particles(i).position);
            particles(i).velocity = inertia + cognitive + social;

            particles(i).position = particles(i).position + particles(i).velocity;
            particles(i).position = max(1, min(numNodes, round(particles(i).position)));

            cost = fitnessFunction(particles(i).position, network, alpha, beta);
            if cost < particles(i).bestCost
                particles(i).bestCost = cost;
                particles(i).bestPosition = particles(i).position;
            end

            if cost < globalBestCost
                globalBestCost = cost;
                globalBestPosition = particles(i).bestPosition;
            end
        end
    end

    Xbest = globalBestPosition;
    bestCost = globalBestCost;
end

% Energy cost calculation function and rounds calculation remain the same as in the original code.
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
    
    cost = alpha *( maxDist/707 )+ beta * maxHop;
end

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




% Simulation Parameters

maxIter = 800;        % Maximum number of iterations
alpha = 0.2;          % Weight for distance
beta = 1 - alpha;     % Weight for hop count

base_station_locations = [
    500, 500;% Scenario 1
    250,250

];

numSensors = 200;  % Number of sensor

% Varying number of sensors from 60 to 140
Gateway_counts = 60:20:140;
round_result = zeros(length(Gateway_counts), size(base_station_locations, 25));
round_results_pso = zeros(length(Gateway_counts), size(base_station_locations, 1));

% Iterate over each base station location
for i = 1:size(base_station_locations, 1)


    for j = 1:length(Gateway_counts)
        numGateways = Gateway_counts(j);
       numParticles = numGateways;    % Number of particles in the swarm
        base_station = base_station_locations(i, :);
        for k = 1 : 7
        % Initialize network and PSO routing
        [network, nextHops, X, positions] = initializeNetwork(numSensors, numGateways, network_range);
        [Xbest, bestCost] = psoRouting(numSensors, numGateways, numParticles, maxIter, alpha, beta, base_station, network_range);
        
        % Calculate energy exhaustion rounds
        rounds = getRoundsUntilFirstNodeExhausted(numSensors, numGateways, base_station, positions, E_elec, epsilon_fs, epsilon_mp, msg_size, d_0, initial_energy, nextHops, network, alpha, beta , Xbest);
        round_result(j,i , k) = rounds;
        end
         round_results_pso(j,i) = mean(round_result(j,i,:));
        % Calculate energy exhausted
       % energy_exhausted = calculateEnergyCost(bestRoute, E_elec, epsilon_fs, epsilon_mp, msg_size, d_0, positions);
       % fprintf('Scenario %d, Sensors %d: Energy exhausted = %.4f J, Rounds = %d\n', i, numSensors, energy_exhausted, rounds);
    end
    end


% Plotting as a bar graph
%{

figure;
bar(sensor_counts, round_results);
xlabel('Number of Sensors');
ylabel('Number of Rounds before Energy Exhaustion');
title('Rounds vs. Sensors for Different Base Station Locations');
legend('Base Station (500, 500)', 'Base Station (250, 0)', 'Base Station (400, 100)', 'Location', 'northwest');
grid on;
%}

% After running pso_final, save the results
save('pso_final_results.mat', 'round_results_pso');


