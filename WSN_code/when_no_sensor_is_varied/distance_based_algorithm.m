clc;
clear all;
rng(42);

% Define the function at the beginning of the script
function [x, y] = Getroundnum(numSensors, numGateways, base_station)
    % Parameters
    
    
    range = 500;
    bandwidth = 1e6;  % 1 MB/s
    msg_size = 525;   % 500 bytes + 25-byte header
    E_elec = 50e-9;  % Energy per bit for electronics
    epsilon_fs = 10e-12;  % Free space energy loss coefficient
    epsilon_mp = 0.0013e-12;  % Multipath energy loss coefficient
    d_0 = 75;
    initial_energy = 10;  % Initial energy level

    % Generate random sensor and gateway coordinates
    rng(42);  % For reproducibility
    sensor_x = rand(1, numSensors) * range;
    sensor_y = rand(1, numSensors) * range;
    gateway_x = rand(1, numGateways) * range;
    gateway_y = rand(1, numGateways) * range;

    % Extract base station coordinates
    base_station_x = base_station(1);
    base_station_y = base_station(2);

    % Calculate distances from gateways to base station
    gateway_distance_to_base_station = sqrt((gateway_x - base_station_x).^2 + (gateway_y - base_station_y).^2);

    % Find the gateway closest to the base station
    [~, base_station_gateway] = min(gateway_distance_to_base_station);

    % Find the nearest gateway for each sensor
    sensor_nearest_gateway = zeros(1, numSensors);
    for i = 1:numSensors
        distances = sqrt((gateway_x - sensor_x(i)).^2 + (gateway_y - sensor_y(i)).^2);
        [~, sensor_nearest_gateway(i)] = min(distances);
    end

    % Create a graph with energy costs
    gateway_energy_cost = calculate_energy_costs(gateway_x, gateway_y, E_elec, epsilon_fs, epsilon_mp, msg_size, d_0);
    G = digraph(gateway_energy_cost);

    % Initialize energy levels for gateways
    gateway_energy_levels = ones(1, numGateways) * initial_energy;

    % Simulate multiple rounds of data transmission
    numRounds = 1000;
    for roundNum = 1:numRounds
        for sensorIndex = 1:numSensors
            start_gateway = sensor_nearest_gateway(sensorIndex);

            % Find the shortest path from the sensor's nearest gateway to the base station gateway
            [path, energy_cost] = shortestpath(G, start_gateway, base_station_gateway);

            % Check if any gateway on the path has insufficient energy
            pathGateways = unique(path);
            insufficient_energy = any(gateway_energy_levels(pathGateways) < energy_cost);

            if insufficient_energy
                x = roundNum;
                y = sensorIndex +1;
                return
            end

            % Deduct energy cost from all gateways along the path
            for g = pathGateways
                gateway_energy_levels(g) = gateway_energy_levels(g) - energy_cost;
            end
        end
    end

    % If no energy is exhausted, set roundNum to the maximum
    x = numRounds;
    y = -1;  % No sensor caused energy exhaustion
end

% Function to calculate energy costs
function energy_costs = calculate_energy_costs(gateway_x, gateway_y, E_elec, epsilon_fs, epsilon_mp, msg_size, d_0)
    numGateways = length(gateway_x);
    energy_costs = inf(numGateways, numGateways);  % Initialize with high values

    for i = 1:numGateways
        for j = 1:numGateways
            if i == j
                continue;
            end
            
            d = sqrt((gateway_x(j) - gateway_x(i))^2 + (gateway_y(j) - gateway_y(i))^2);
            
            if d < d_0
                energy_transmitter = E_elec * msg_size + epsilon_fs * (d^2) * msg_size;
            else
                energy_transmitter = E_elec * msg_size + epsilon_mp * (d^4) * msg_size;
            end
            
            energy_receiver = E_elec * msg_size;
            
            energy_costs(i, j) = energy_transmitter + energy_receiver;
        end
    end
end

% Now you can call the function after defining it

base_station_locations = [
    500, 500;
    250,250];  % Scenario 1
 

numGateways = 60;  % Number of sensors

% Vary the number of gateways and collect data for each scenario
sensor_counts = 100:100:700;  % Varying number of gateways from 10 to 50
round_results = zeros(length(sensor_counts), size(base_station_locations, 1));

for i = 1:size(base_station_locations, 1)
    for j = 1:length(sensor_counts)  % For each gateway count
        numSensors = sensor_counts(j);
        
        base_station = base_station_locations(i, :);
        % Simulate and get the number of rounds before energy is exhausted
        [round_results_dist(j,i), ~] = Getroundnum(numSensors, numGateways, base_station);
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
save('distance_based_algorithm_final.mat', 'round_results_dist');