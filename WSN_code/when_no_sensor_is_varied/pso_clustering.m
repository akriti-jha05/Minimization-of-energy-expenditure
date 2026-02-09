clc;
clear all;
close all;

% Network Parameters
network_range = 500;      % Range of the network
numClusters = 4;          % Number of clusters
msg_size = 525;           % Message size (500 bytes + 25-byte header)
E_elec = 50e-9;           % Energy per bit for electronics
epsilon_fs = 10e-12;      % Free space energy loss coefficient
epsilon_mp = 0.0013e-12;  % Multipath energy loss coefficient
d_0 = 75;                 % Threshold distance for free space/multipath
initial_energy = 10;       % Initial energy for each node (J)

% Simulation Parameters
numParticles = 10;        % Number of particles for PSO
maxIter = 100;            % Maximum iterations for PSO
alpha = 0.5;              % Weight for communication distance
beta = 1 - alpha;         % Weight for residual energy

% Define Base Station Scenarios
base_station_locations = [
    250, 250;  % Scenario 1: Sink in the middle of the sensing field
    500, 500   % Scenario 2: Sink at the top-center of the sensing field
];

% Number of Gateways and Sensor Nodes
numGateways = 60;
sensor_counts = 100:100:700; % Sensor nodes from 100 to 700

% Preallocate results
round_results_pso = zeros(length(sensor_counts), size(base_station_locations, 1));

% Main Simulation Loop
for scenario = 1:size(base_station_locations, 1)
    base_station = base_station_locations(scenario, :);

    for sensorIdx = 1:length(sensor_counts)
        numSensors = sensor_counts(sensorIdx);

        % Initialize Network and Nodes
        [network, positions] = initializeNetwork(numSensors, numGateways, network_range);

        % Run PSO to determine optimal cluster heads
        [clusterHeads, bestCost] = psoClustering(numSensors, numGateways, numParticles, maxIter, network, positions, base_station, alpha, beta, numClusters);

        % Calculate rounds until the first node exhausts its energy
        rounds = simulateRoundsUntilDepletion(numSensors, numGateways, clusterHeads, positions, base_station, E_elec, epsilon_fs, epsilon_mp, msg_size, d_0, initial_energy);
        round_results_pso(sensorIdx, scenario) = rounds;
    end
end

% Plot Results
figure;
bar(sensor_counts, round_results_pso);
xlabel('Number of Sensors');
ylabel('Rounds until First Node Exhaustion');
title('Performance of PSO-based Clustering');
legend('Sink at Center', 'Sink at Top Center', 'Location', 'northwest');
grid on;

% Save Results
save('pso_clustering_results.mat', 'round_results_pso');

%% Functions

% Initialize the network and node positions
function [network, positions] = initializeNetwork(numSensors, numGateways, network_range)
    totalNodes = numSensors + numGateways;
    positions = rand(totalNodes, 2) * network_range; % Random positions in the field
    network = zeros(totalNodes, totalNodes);

    % Define communication range (120 meters for sensor-to-gateway, gateway-to-gateway)
    comm_range = 120;

    % Fill adjacency matrix for connectivity
    for i = 1:totalNodes
        for j = 1:totalNodes
            if i ~= j
                d = sqrt((positions(i, 1) - positions(j, 1))^2 + (positions(i, 2) - positions(j, 2))^2);
                if d <= comm_range
                    network(i, j) = d;
                end
            end
        end
    end
end

% Particle Swarm Optimization for clustering
function [clusterHeads, bestCost] = psoClustering(numSensors, numGateways, numParticles, maxIter, network, positions, baseStation, alpha, beta, numClusters)
    totalNodes = numSensors + numGateways;

    % Initialize particles
    particles = struct();
    for i = 1:numParticles
        particles(i).position = randi([1, totalNodes], 1, numClusters); % Random cluster heads
        particles(i).velocity = zeros(1, numClusters); % Initial velocity
        particles(i).bestPosition = particles(i).position; % Initial best position
        particles(i).bestCost = fitnessFunction(particles(i).position, network, positions, baseStation, alpha, beta); % Initial cost
    end

    % Initialize global best
    [globalBestCost, idx] = min([particles.bestCost]);
    globalBestPosition = particles(idx).bestPosition;

    % PSO parameters
    w = 0.5;  % Inertia weight
    c1 = 1.5; % Cognitive component
    c2 = 1.5; % Social component

    % PSO main loop
    for iter = 1:maxIter
        for i = 1:numParticles
            % Update velocity
            inertia = w * particles(i).velocity;
            cognitive = c1 * rand * (particles(i).bestPosition - particles(i).position);
            social = c2 * rand * (globalBestPosition - particles(i).position);
            particles(i).velocity = inertia + cognitive + social;

            % Update position
            particles(i).position = round(particles(i).position + particles(i).velocity);
            particles(i).position = max(1, min(totalNodes, particles(i).position)); % Bound positions

            % Evaluate fitness
            cost = fitnessFunction(particles(i).position, network, positions, baseStation, alpha, beta);
            if cost < particles(i).bestCost
                particles(i).bestCost = cost;
                particles(i).bestPosition = particles(i).position;
            end

            % Update global best
            if cost < globalBestCost
                globalBestCost = cost;
                globalBestPosition = particles(i).position;
            end
        end
    end

    clusterHeads = globalBestPosition;
    bestCost = globalBestCost;
end

% Fitness function for PSO
function cost = fitnessFunction(clusterHeads, network, positions, baseStation, alpha, beta)
    % Communication distance and energy cost
    totalDistance = 0;
    for i = 1:length(clusterHeads)
        totalDistance = totalDistance + sqrt((positions(clusterHeads(i), 1) - baseStation(1))^2 + (positions(clusterHeads(i), 2) - baseStation(2))^2);
    end

    % Residual energy (assumed to be uniform initially)
    residualEnergy = 1; % Placeholder for uniform initial energy

    % Cost as weighted sum
    cost = alpha * totalDistance + beta * residualEnergy;
end

% Simulate rounds until the first node exhausts energy
function rounds = simulateRoundsUntilDepletion(numSensors, numGateways, clusterHeads, positions, baseStation, E_elec, epsilon_fs, epsilon_mp, msg_size, d_0, initial_energy)
    totalNodes = numSensors + numGateways;
    energyLevels = initial_energy * ones(1, totalNodes); % Initialize energy levels
    rounds = 0;

    while all(energyLevels > 0)
        for i = 1:numSensors
            % Select nearest cluster head
            distances = arrayfun(@(ch) sqrt((positions(i, 1) - positions(ch, 1))^2 + (positions(i, 2) - positions(ch, 2))^2), clusterHeads);
            [~, idx] = min(distances);
            selectedCH = clusterHeads(idx);

            % Energy cost for transmission
            d = distances(idx);
            if d < d_0
                E_tx = E_elec * msg_size + epsilon_fs * msg_size * d^2;
            else
                E_tx = E_elec * msg_size + epsilon_mp * msg_size * d^4;
            end

            % Update energy levels
            energyLevels(i) = energyLevels(i) - E_tx;
            energyLevels(selectedCH) = energyLevels(selectedCH) - (E_elec * msg_size);

            % Check if any node is depleted
            if energyLevels(i) <= 0 || energyLevels(selectedCH) <= 0
                return;
            end
        end
        rounds = rounds + 1;
    end
end
