clc;
clear all;
close all;

% Network Parameters
network_range = 500;      % Range of the network
numSensors = 200;         % Fixed number of sensors
msg_size = 525;           % Message size (500 bytes + 25-byte header)
E_elec = 50e-9;           % Energy per bit for electronics
epsilon_fs = 10e-12;      % Free space energy loss coefficient
epsilon_mp = 0.0013e-12;  % Multipath energy loss coefficient
d_0 = 75;                 % Threshold distance for free space/multipath
initial_energy = 2;       % Initial energy for each node (J)

% Simulation Parameters
numParticles = 10;        % Number of particles for PSO
maxIter = 200;             % Reduced max iterations for PSO
alpha = 0.5;              % Weight for communication distance
beta = 1 - alpha;         % Weight for residual energy

% Base Station Locations
base_station_locations = [
    250, 250;  % Base Station 1
    500, 500   % Base Station 2
];

% Varying Number of Gateways
gateway_counts = 60:20:200; % Gateways from 60 to 200 in steps of 20
round_results_pso = zeros(length(gateway_counts), size(base_station_locations, 1)); % Preallocate results

% Main Simulation Loop
for baseStationIdx = 1:size(base_station_locations, 1)
    base_station = base_station_locations(baseStationIdx, :);

    for gatewayIdx = 1:length(gateway_counts)
        numGateways = gateway_counts(gatewayIdx);

        % Initialize Network and Nodes
        [network, positions, distanceMatrix] = initializeNetwork(numSensors, numGateways, network_range);

        % Run PSO to determine optimal cluster heads
        [clusterHeads, bestCost] = psoClustering(numSensors, numGateways, numParticles, maxIter, distanceMatrix, positions, base_station, alpha, beta, numGateways);

        % Calculate rounds until the first node exhausts its energy
        rounds = simulateRoundsUntilDepletion(numSensors, numGateways, clusterHeads, positions, base_station, E_elec, epsilon_fs, epsilon_mp, msg_size, d_0, initial_energy, distanceMatrix);
        round_results_pso(gatewayIdx, baseStationIdx) = rounds;
    end
end

% Plot Results
figure;
bar(gateway_counts, round_results_pso, 'grouped');
xlabel('Number of Gateways');
ylabel('Rounds until First Node Exhaustion');
title('Performance of PSO-based Clustering with Varying Gateways');
legend('Base Station (250, 250)', 'Base Station (500, 500)', 'Location', 'northwest');
grid on;

% Save Results
save('pso_optimized_vary_gateways_results.mat', 'round_results_pso');

%% Functions

% Initialize the network and node positions
function [network, positions, distanceMatrix] = initializeNetwork(numSensors, numGateways, network_range)
    totalNodes = numSensors + numGateways;
    positions = rand(totalNodes, 2) * network_range; % Random positions in the field
    network = zeros(totalNodes, totalNodes);

    % Precompute distance matrix
    distanceMatrix = zeros(totalNodes);
    for i = 1:totalNodes
        for j = i+1:totalNodes
            d = sqrt((positions(i, 1) - positions(j, 1))^2 + (positions(i, 2) - positions(j, 2))^2);
            distanceMatrix(i, j) = d;
            distanceMatrix(j, i) = d; % Symmetric
        end
    end

    % Define communication range (120 meters for connectivity)
    comm_range = 100;
    network(distanceMatrix <= comm_range) = 1; % Mark connections
end

% Particle Swarm Optimization for clustering
function [clusterHeads, bestCost] = psoClustering(numSensors, numGateways, numParticles, maxIter, distanceMatrix, positions, baseStation, alpha, beta, numClusters)
    totalNodes = numSensors + numGateways;

    % Initialize particles
    particles = struct();
    for i = 1:numParticles
        particles(i).position = randi([1, totalNodes], 1, numClusters); % Random cluster heads
        particles(i).velocity = zeros(1, numClusters); % Initial velocity
        particles(i).bestPosition = particles(i).position; % Initial best position
        particles(i).bestCost = fitnessFunction(particles(i).position, distanceMatrix, positions, baseStation, alpha, beta); % Initial cost
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
            cost = fitnessFunction(particles(i).position, distanceMatrix, positions, baseStation, alpha, beta);
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
function cost = fitnessFunction(clusterHeads, distanceMatrix, positions, baseStation, alpha, beta)
    % Communication distance and energy cost
    totalDistance = 0;
    for i = 1:length(clusterHeads)
        totalDistance = totalDistance + sqrt((positions(clusterHeads(i), 1) - baseStation(1))^2 + (positions(clusterHeads(i), 2) - baseStation(2))^2);
    end

    % Residual energy (placeholder for uniform initial energy)
    residualEnergy = 1; % Uniform

    % Cost as weighted sum
    cost = alpha * totalDistance + beta * residualEnergy;
end

% Simulate rounds until the first node exhausts energy
function rounds = simulateRoundsUntilDepletion(numSensors, numGateways, clusterHeads, positions, baseStation, E_elec, epsilon_fs, epsilon_mp, msg_size, d_0, initial_energy, distanceMatrix)
    totalNodes = numSensors + numGateways;
    energyLevels = initial_energy * ones(1, totalNodes); % Initialize energy levels
    rounds = 0;

    while all(energyLevels > 0)
        for i = 1:numSensors
            % Select nearest cluster head
            [~, idx] = min(distanceMatrix(i, clusterHeads));
            selectedCH = clusterHeads(idx);

            % Energy cost for transmission
            d = distanceMatrix(i, selectedCH);
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
