%% Nathan Marchant May 2025
% Written for Pavlovian conflict task
% --- MODIFIED: Re-introduced flexible session grouping for early/late comparisons ---
% --- MODIFIED: Reverted to explicit plotting blocks to fix legend errors ---

% --- Initialization ---
clear all;
close all; % Good practice to close figures

% --- MODIFIED: Define session groups ---
% List the session names to pool for the "early" and "late" groups.
early_sessions = {'shock_1', 'shock_1'};
late_sessions = {'shock_4', 'shock_4'};

% --- Define the folder containing your extracted .mat files ---
dataFolder = 'C:\Photometry\PavConf\DrPhotom_Extracted\Alcohol conditioning V4 photometry_SHOCK'; % <--- SET YOUR FOLDER PATH HERE
files = dir(fullfile(dataFolder, '*.mat'));
fprintf('Found %d data files in: %s\n', length(files), dataFolder);

% --- MODIFIED: Initialize data matrices for each condition and session group ---
CSP_no_shock_early = []; CSP_shock_early = []; CSM_early = [];
CSP_no_shock_late = [];  CSP_shock_late = [];  CSM_late = [];

% --- Data Extraction and Processing Loop ---
for k = 1:length(files)
    currentFile = fullfile(dataFolder, files(k).name);
    fprintf('Loading and processing: %s\n', files(k).name);
    
    load(currentFile); % This loads the 'sesdat' structure
    
    % Check if the loaded file contains the necessary data structure
    if exist('sesdat', 'var') && isfield(sesdat, 'traces_z') && ~isempty(sesdat.traces_z) && isfield(sesdat, 'session')
        
        % Ensure column 3 exists for shock checking
        if size(sesdat.traces_z, 2) < 3
            sesdat.traces_z(:,3) = 0; % Assume no shock if column is missing
        end
        
        % Find indices for each of the three conditions
        idx_csp_no_shock = (sesdat.traces_z(:, 2) == 1 & sesdat.traces_z(:, 3) == 0);
        idx_csp_shock = (sesdat.traces_z(:, 2) == 1 & sesdat.traces_z(:, 3) == 1);
        idx_csm = (sesdat.traces_z(:, 2) == 2);
        
        % Append data to the correct matrix based on session
        if ismember(sesdat.session, early_sessions)
            if any(idx_csp_no_shock), CSP_no_shock_early = [CSP_no_shock_early; sesdat.traces_z(idx_csp_no_shock, 5:end)]; end
            if any(idx_csp_shock),    CSP_shock_early    = [CSP_shock_early;    sesdat.traces_z(idx_csp_shock, 5:end)]; end
            if any(idx_csm),          CSM_early          = [CSM_early;          sesdat.traces_z(idx_csm, 5:end)]; end
        elseif ismember(sesdat.session, late_sessions)
            if any(idx_csp_no_shock), CSP_no_shock_late = [CSP_no_shock_late; sesdat.traces_z(idx_csp_no_shock, 5:end)]; end
            if any(idx_csp_shock),    CSP_shock_late    = [CSP_shock_late;    sesdat.traces_z(idx_csp_shock, 5:end)]; end
            if any(idx_csm),          CSM_late          = [CSM_late;          sesdat.traces_z(idx_csm, 5:end)]; end
        end
        
    else
        fprintf('  -> Skipping file. "sesdat.traces_z" or "session" field is missing or empty.\n');
    end
    clear sesdat;
end

fprintf('\nData collation complete.\n');

% --- Data Cleaning for each group ---
% Function to clean data (removes outliers and NaNs)
clean_data = @(data) data(~(any(isnan(data), 2) | abs(mean(data, 2)) > 20), :);

CSP_no_shock_early = clean_data(CSP_no_shock_early);
CSP_shock_early = clean_data(CSP_shock_early);
CSM_early = clean_data(CSM_early);
CSP_no_shock_late = clean_data(CSP_no_shock_late);
CSP_shock_late = clean_data(CSP_shock_late);
CSM_late = clean_data(CSM_late);

fprintf('\nData cleaning complete.\n');

% --- Plotting Section ---
% Define time vector (find first non-empty matrix to base it on)
all_data = {CSP_no_shock_early, CSP_shock_early, CSM_early, CSP_no_shock_late, CSP_shock_late, CSM_late};
time = [];
for i = 1:length(all_data)
    if ~isempty(all_data{i})
        time = linspace(-10, 20, size(all_data{i}, 2));
        break;
    end
end
if isempty(time), error('No data available for plotting.'); end

% Define colors
green = [0.10, 0.60, 0.00]; red = [1.0, 0, 0]; orange = [1.00, 0.45, 0.00];
light_green = [0.5, 0.8, 0.5]; light_red = [1.0, 0.5, 0.5]; light_orange = [1.0, 0.7, 0.4];
black1 = [0.75,0.75,0.75]; black2 = [0.5,0.5,0.5]; black3 = [0.25,0.25,0.25];

% --- Statistical Parameters ---
p = 0.01;
thres = 8;

%% ----------------- Plot 1: no shock sessions comparison -----------------
% --- Plot 1: no shock Early vs Late ---
if ~isempty(CSP_no_shock_early) && ~isempty(CSP_no_shock_late)
    plot_title = 'CS no shock: Early vs Late';
    labels = {sprintf('CSP (n=%d)', size(CSP_no_shock_early,1)), sprintf('CSM (n=%d)', size(CSP_no_shock_late,1))};
    create_comparison_plot(CSP_no_shock_early, CSP_no_shock_late, labels, {green, orange}, plot_title, time);
end

%% ----------------- Plot 2: shock sessions comparison -----------------
% --- Plot 2: Shock Early vs Late ---
if ~isempty(CSP_shock_early) && ~isempty(CSP_shock_late)
    plot_title = 'CS Shock: Early vs Late';
    labels = {sprintf('CSP (n=%d)', size(CSP_shock_early,1)), sprintf('CSM (n=%d)', size(CSP_shock_late,1))};
    create_comparison_plot(CSP_shock_early, CSP_shock_late, labels, {green, orange}, plot_title, time);
end

%% ----------------- Plot 3: csM sessions comparison -----------------
% --- Plot 3: CSM  Early vs Late ----
if ~isempty(CSM_early) && ~isempty(CSM_late)
    plot_title = 'CSM: Early vs Late';
    labels = {sprintf('CSP (n=%d)', size(CSM_early,1)), sprintf('CSM (n=%d)', size(CSM_late,1))};
    create_comparison_plot(CSM_early, CSM_late, labels, {green, orange}, plot_title, time);
end

%% ----------------- Plot 4: Early session -----------------
% --- Plot 4: Early ----
if ~isempty(CSP_no_shock_early) && ~isempty(CSP_shock_early) && ~isempty(CSM_early)
    plot_title = 'Early';
    labels = {sprintf('CSP (n=%d)', size(CSP_no_shock_early,1)), sprintf('CSM (n=%d)', size(CSP_shock_early,1)), sprintf('CSM (n=%d)', size(CSM_early,1))};
    create_three_comparison_plot(CSP_no_shock_early, CSP_shock_early, CSM_early, labels, {green, orange, red}, plot_title, time);
end

%% ----------------- Plot 5: Latesession -----------------

if ~isempty(CSP_no_shock_late) && ~isempty(CSP_shock_late) && ~isempty(CSM_late)
    plot_title = 'Late';
    labels = {sprintf('CSP (n=%d)', size(CSP_no_shock_late,1)), sprintf('CSM (n=%d)', size(CSP_shock_late,1)), sprintf('CSM (n=%d)', size(CSM_late,1))};
    create_three_comparison_plot(CSP_no_shock_late, CSP_shock_late, CSM_late, labels, {green, orange, red}, plot_title, time);
end

%% ----------------- Plotting Function -----------------
function create_comparison_plot(data1, data2, labels, colors, plot_title, time)
    % Creates a comparison plot with two traces, shaded error, and stats.
    
    figure('Name', plot_title, 'NumberTitle', 'off');
    hold on;
    
    % --- Stats ---
    [perm_test, ~] = permTest_array(data1, data2, 1000);
    
    % --- Plotting ---
    datasets = {data1, data2};
    for i = 1:2
        ds = datasets{i};
        mean_trace = mean(ds, 1);
        sem_trace = std(ds, 0, 1) / sqrt(size(ds, 1));
        
        plot(time, mean_trace, 'Color', colors{i}, 'LineWidth', 1.5);
        jbfill(time, mean_trace - sem_trace, mean_trace + sem_trace, colors{i}, 'none', 0, 0.2);
    end
    
    % --- Set Y-axis limits and position significance bar ---
    ylim([-2 12]); % Hard-code the axis range
    y_lim = ylim; % Get the new limits
    
    p_val = 0.01;
    thres = 8;
    % Position the significance bar just inside the bottom of the new fixed axis
    sig_y = y_lim(1) + 0.5; 
    
    tmp = find(perm_test(1, :) < p_val);
    id = tmp(consec_idx(tmp, thres));
    if ~isempty(id)
        plot(time(id), sig_y * ones(size(time(id))), 's', 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'Color', 'k');
    end
    
    % --- Final Touches ---
    legend(labels, 'Location', 'northwest');
    title(plot_title);
    xlabel('Time from Cue Onset (s)');
    ylabel('Z-Score');
    line([0 0], y_lim, 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
    line([10 10], y_lim, 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
    line(xlim, [0 0], 'Color', 'k', 'LineStyle', '-');
    box on;
    hold off;

end

%% ----------------- Plotting Function -----------------
function create_three_comparison_plot(data1, data2, data3, labels, colors, plot_title, time)
    % Creates a comparison plot with THREE traces, shaded error, and stats.
    
    figure('Name', plot_title, 'NumberTitle', 'off');
    hold on;
    
    % --- Stats for all three pairs ---
    [perm_test_1v2, ~] = permTest_array(data1, data2, 1000);
    [perm_test_1v3, ~] = permTest_array(data1, data3, 1000);
    [perm_test_2v3, ~] = permTest_array(data2, data3, 1000);
    
    % --- Plotting ---
    datasets = {data1, data2, data3};
    for i = 1:3
        ds = datasets{i};
        if isempty(ds), continue; end % Skip if data is empty
        
        mean_trace = mean(ds, 1);
        sem_trace = std(ds, 0, 1) / sqrt(size(ds, 1));
        
        plot(time, mean_trace, 'Color', colors{i}, 'LineWidth', 1.5);
        jbfill(time, mean_trace - sem_trace, mean_trace + sem_trace, colors{i}, 'none', 0, 0.2);
    end
    
    % --- Set Y-axis limits and position significance bars ---
    ylim([-2 12]); % Hard-code the axis range
    y_lim = ylim; % Get the new limits
    
    p_val = 0.01;
    thres = 8;
    
    % Position the significance bars at different heights
    sig_y_1 = y_lim(1) + 0.5; 
    sig_y_2 = y_lim(1) + 0.65;
    sig_y_3 = y_lim(1) + 0.8;
    
    % Plot significance for 1 vs 2
    tmp = find(perm_test_1v2(1, :) < p_val);
    id = tmp(consec_idx(tmp, thres));
    if ~isempty(id)
        plot(time(id), sig_y_1 * ones(size(time(id))), 's', 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'Color', 'k');
    end
    
    % Plot significance for 1 vs 3
    tmp = find(perm_test_1v3(1, :) < p_val);
    id = tmp(consec_idx(tmp, thres));
    if ~isempty(id)
        plot(time(id), sig_y_2 * ones(size(time(id))), 's', 'MarkerSize', 7, 'MarkerFaceColor', [0.4 0.4 0.4], 'Color', [0.4 0.4 0.4]);
    end
    
    % Plot significance for 2 vs 3
    tmp = find(perm_test_2v3(1, :) < p_val);
    id = tmp(consec_idx(tmp, thres));
    if ~isempty(id)
        plot(time(id), sig_y_3 * ones(size(time(id))), 's', 'MarkerSize', 7, 'MarkerFaceColor', [0.7 0.7 0.7], 'Color', [0.7 0.7 0.7]);
    end
    
    % --- Final Touches ---
    legend(labels, 'Location', 'northwest');
    title(plot_title);
    xlabel('Time from Cue Onset (s)');
    ylabel('Z-Score');
    line([0 0], y_lim, 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
    line([10 10], y_lim, 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
    line(xlim, [0 0], 'Color', 'k', 'LineStyle', '-');
    box on;
    hold off;
end
