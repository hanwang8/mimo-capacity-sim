% Capacity of MIMO Rayleigh fading channels
% Simulation: Ergodic and 99% Outage Capacity for various antenna configurations
clc;
clear;
close all;
rng(1);                       % Reproducible results

%% Simulation Parameters
SNR_dB = 0:2:40;             % SNR range in dB
SNR_lin = 10.^(SNR_dB/10);   % Convert to linear scale
num_snr_points = length(SNR_dB);

N_iter = 1e4;                % Monte Carlo iterations (10,000)
P_out_target = 0.01;         % Target outage probability (1% for 99% capacity)

% Antenna configurations [Tx, Rx]
antenna_configs = [1 1; 1 2; 2 1; 2 2];
num_configs = size(antenna_configs, 1);

% Pre-allocate storage for results
ergodic_capacity = zeros(num_snr_points, num_configs);
outage_capacity = zeros(num_snr_points, num_configs);
legend_str = cell(num_configs, 1);

disp('--- Starting MIMO Capacity Simulation ---');

%% Main Simulation Loop
for c_idx = 1:num_configs
    M = antenna_configs(c_idx, 1); % Tx antennas
    N = antenna_configs(c_idx, 2); % Rx antennas
    
    % Store legend entry
    legend_str{c_idx} = sprintf('M=%d, N=%d', M, N);
    fprintf('Simulating Configuration: %d x %d ...\n', M, N);
    
    % Optimization: Pre-define identity matrices outside the loop
    I_M = eye(M);
    I_N = eye(N);

    for s_idx = 1:num_snr_points
        rho = SNR_lin(s_idx); % Current SNR
        C_inst = zeros(N_iter, 1); % Store instantaneous capacities
        
        for k = 1:N_iter
            % Generate Channel Matrix H ~ CN(0,1)
            % Normalization by sqrt(2) ensures unit variance per entry
            H = (randn(N,M) + 1j*randn(N,M)) / sqrt(2);
            
            % Calculate Capacity
            % Optimization: Use Sylvester's determinant identity to compute
            % determinant of the smaller dimension matrix.
            if M <= N
                W = H' * H;
                % C = log2(det(I + (rho/M) * H'*H))
                C_inst(k) = real(log2(det(I_M + (rho/M) * W)));
            else
                W = H * H';
                % C = log2(det(I + (rho/M) * HH'))
                C_inst(k) = real(log2(det(I_N + (rho/M) * W)));
            end
        end
        
        % 1. Calculate Ergodic Capacity (Mean)
        ergodic_capacity(s_idx, c_idx) = mean(C_inst);
        
        % 2. Calculate 99% Outage Capacity
        % The capacity value that 99% of channels can support
        sorted_caps = sort(C_inst);
        idx_out = ceil(P_out_target * N_iter); 
        % Note: If P_out is 1%, we look at the 1st percentile of sorted data
        outage_capacity(s_idx, c_idx) = sorted_caps(idx_out);
    end
end

disp('--- Simulation Finished ---');

%% Plotting Results

% Figure 1: 99% Outage Capacity
fig1 = figure;
plot(SNR_dB, outage_capacity, '-s', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('SNR (dB)');
ylabel('Capacity (bits/s/Hz)');
title('99% Outage Capacity of MIMO Fading Channels');
legend(legend_str, 'Location', 'NorthWest');

% Figure 2: Ergodic Capacity
fig2 = figure;
plot(SNR_dB, ergodic_capacity, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('SNR (dB)');
ylabel('Capacity (bits/s/Hz)');
title('Ergodic Capacity of MIMO Fading Channels');
legend(legend_str, 'Location', 'NorthWest');

%% Save the results
resultsDir = 'results';
if ~exist(resultsDir, 'dir'); mkdir(resultsDir); end
labels = {'1x1', '1x2', '2x1', '2x2'};
names = [{'SNR_dB'}, strcat('ergodic_', labels), strcat('outage1pct_', labels)];
results = array2table([SNR_dB(:), ergodic_capacity, outage_capacity], 'VariableNames', names);
writetable(results, fullfile(resultsDir, 'capacity_results.csv'));
exportgraphics(fig1, fullfile(resultsDir, 'outage_capacity.png'), 'Resolution', 150);
exportgraphics(fig2, fullfile(resultsDir, 'ergodic_capacity.png'), 'Resolution', 150);
