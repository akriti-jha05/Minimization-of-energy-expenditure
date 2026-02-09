clc;
close all;
clear all;

% Load data from each algorithm
data1 = load('distance_based_algorithm_final.mat'); % Ensure this file contains 'round_results_dist'
data2 = load('pso_final_results.mat'); % Ensure this file contains 'round_results_pso'
data3 = load('grey_wolf.mat'); % Ensure this file contains 'round_results_gwo'

data5 = load('cray_optimizer.mat'); % Ensure this file contains 'round_results_cray_optimizer'

% Extract the results from the loaded data
round_results_dist = data1.round_results_dist; 
round_results_pso = data2.round_results_pso;
round_results_gwo = data3.round_results_gwo;
round_results_cray_optimizer = data5.round_results_cray_optimizer;

% Number of gateways varying from 60 to 140
Gateways_counts = 60:20:140;

% Combine the results into matrices for each base station location
round_results_bs1 = [round_results_dist(:,1), round_results_pso(:,1),  round_results_gwo(:,1),round_results_cray_optimizer(:,1) ];
round_results_bs2 = [round_results_dist(:,2), round_results_pso(:,2), round_results_gwo(:,2),round_results_cray_optimizer(:,2) ];


% Plot for base station at [500, 500]
figure;
bar(Gateways_counts, round_results_bs1, 'grouped'); % Since round_results_bs1 is already a 5xN matrix, no transpose is needed
grid on;
xlabel('Number of Gateways');
ylabel('Number of Rounds before first node die');
title('Comparison of Algorithms (Base Station at [500,500])');
legend('Distance-based', 'PSO',  'GW0' , 'COA');

% Plot for base station at [250, 250]
figure;
bar(Gateways_counts, round_results_bs2, 'grouped'); % Since round_results_bs2 is already a 5xN matrix, no transpose is needed
grid on;
xlabel('Number of Gateways');
ylabel('Number of Rounds before First Node die');
title('Comparison of Algorithms (Base Station at [250,250])');
legend('Distance-based', 'PSO',  'GWO' ,'COA');
