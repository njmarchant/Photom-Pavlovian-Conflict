%% Nathan Marchant May 2025
% Written for Pavlovian conflict task
% Includes: flexible grouping for first and final sessions ---
% Includes: Integrated bootstrapping tests into plotting function ---
% Includes: Stored permutation tests in a named struct ---
% Includes: output for significance windows ---

% --- Initialization ---
clear all;
close all; % Good practice to close figures

%% ---  Define which sessions to group together ---
% List the session names you want to pool for the "first" group
first_sessions = {'retrain_1', 'retrain_2'}; 
% List the session names for the "final" group for most rats
final_sessions_default = {'retrain_4', 'retrain_5'}; 
% List the session names for the "final" group for special case rats
final_sessions_special = {'retrain_4'}; 
special_rats = {'R15', 'R16'}; % List of rats that use the special final group

% --- Define the folder containing your extracted .mat files ---
dataFolder = 'C:\Photometry\PavConf\DrPhotom_Extracted\Alcohol conditioning V4 photometry'; % <--- SET YOUR FOLDER PATH HERE
files = dir(fullfile(dataFolder, '*.mat'));
fprintf('Found %d data files in: %s\n', length(files), dataFolder);

% --- Initialize data matrices for first and final sessions ---
CSP_first = []; CSM_first = [];
CSP_final = []; CSM_final = [];

% --- Data Extraction and Processing Loop ---
for k = 1:length(files)
    currentFile = fullfile(dataFolder, files(k).name);
    fprintf('Loading and processing: %s\n', files(k).name);
    
    rat_token = regexp(currentFile, 'data_(R\d+)\.mat$', 'tokens', 'once');
    final_sessions_to_use = final_sessions_default;
    if ~isempty(rat_token) && ismember(rat_token{1}, special_rats)
        final_sessions_to_use = final_sessions_special;
    end
    
    load(currentFile);
    
    if exist('sesdat', 'var') && isfield(sesdat, 'traces_z') && ~isempty(sesdat.traces_z) && isfield(sesdat, 'session')
        csp_indices = sesdat.traces_z(:, 2) == 1;
        csm_indices = sesdat.traces_z(:, 2) == 2;
        
        if ismember(sesdat.session, first_sessions)
            if any(csp_indices), CSP_first = [CSP_first; sesdat.traces_z(csp_indices, 5:end)]; end
            if any(csm_indices), CSM_first = [CSM_first; sesdat.traces_z(csm_indices, 5:end)]; end
        elseif ismember(sesdat.session, final_sessions_to_use)
            if any(csp_indices), CSP_final = [CSP_final; sesdat.traces_z(csp_indices, 5:end)]; end
            if any(csm_indices), CSM_final = [CSM_final; sesdat.traces_z(csm_indices, 5:end)]; end
        end
    else
        fprintf('  -> Skipping file. "sesdat" variable, "traces_z", or "session" field is missing or empty.\n');
    end
    clear sesdat;
end
fprintf('\nData collation complete.\n');

% --- Data Cleaning for each group ---
clean_data = @(data) data(~(any(isnan(data), 2) | abs(mean(data, 2)) > 20), :);
CSP_first = clean_data(CSP_first);
CSM_first = clean_data(CSM_first);
CSP_final = clean_data(CSP_final);
CSM_final = clean_data(CSM_final);
fprintf('\nData cleaning complete.\n');
fprintf('Total First Sessions CSP trials: %d\n', size(CSP_first, 1));
fprintf('Total First Sessions CSM trials: %d\n', size(CSM_first, 1));
fprintf('Total Final Sessions CSP trials: %d\n', size(CSP_final, 1));
fprintf('Total Final Sessions CSM trials: %d\n', size(CSM_final, 1));

% --- Perform Bootstrapping on all datasets ---
fprintf('\nPerforming bootstrapping tests...\n');
btsrp = struct();
% NOTE: Assumes 'bootstrap_data' and 'CIadjust' are in your path.
tmp = bootstrap_data(CSP_first, 5000, 0.001);
btsrp.csp_first = CIadjust(tmp(1,:),tmp(2,:),tmp,size(CSP_first, 1),2);
tmp = bootstrap_data(CSM_first, 5000, 0.001);
btsrp.csm_first = CIadjust(tmp(1,:),tmp(2,:),tmp,size(CSM_first, 1),2);
tmp = bootstrap_data(CSP_final, 5000, 0.001);
btsrp.csp_final = CIadjust(tmp(1,:),tmp(2,:),tmp,size(CSP_final, 1),2);
tmp = bootstrap_data(CSM_final, 5000, 0.001);
btsrp.csm_final = CIadjust(tmp(1,:),tmp(2,:),tmp,size(CSM_final, 1),2);
fprintf('Bootstrapping complete.\n');

% --- Perform Permutation Tests and store in a struct ---
fprintf('\nPerforming permutation tests...\n');
perm = struct();
% NOTE: Assumes 'permTest_array' is in your path.
if ~isempty(CSP_first) && ~isempty(CSM_first), [perm.first_csp_vs_csm, ~] = permTest_array(CSP_first, CSM_first, 1000); end
if ~isempty(CSP_final) && ~isempty(CSM_final), [perm.final_csp_vs_csm, ~] = permTest_array(CSP_final, CSM_final, 1000); end
if ~isempty(CSP_first) && ~isempty(CSP_final), [perm.csp_first_vs_final, ~] = permTest_array(CSP_first, CSP_final, 1000); end
if ~isempty(CSM_first) && ~isempty(CSM_final), [perm.csm_first_vs_final, ~] = permTest_array(CSM_first, CSM_final, 1000); end
fprintf('Permutation tests complete.\n');

% --- Define Time Vector and Colors ---
if isempty(CSP_first) && isempty(CSM_first) && isempty(CSP_final) && isempty(CSM_final)
    error('No data available for plotting. Check source files.');
end
% Define time vector
if ~isempty(CSP_first), time = linspace(-10, 20, size(CSP_first, 2));
elseif ~isempty(CSM_first), time = linspace(-10, 20, size(CSM_first, 2));
elseif ~isempty(CSP_final), time = linspace(-10, 20, size(CSP_final, 2));
else, time = linspace(-10, 20, size(CSM_final, 2));
end
% Define colors
green = [0.10, 0.60, 0.00]; light_green = [0.5, 0.8, 0.5];
orange = [1.00, 0.45, 0.00]; light_orange = [1.0, 0.7, 0.4];

% --- NEW: Define Significance Parameters ---
p_val = 0.01; % p-value for permutation test
thres = 8;    % Consecutive data points threshold for significance

% --- NEW: Calculate and Store Significance Windows ---
fprintf('\nCalculating significance windows...\n');
sig_windows = struct();
% NOTE: Assumes 'consec_idx' is in your path.

% Permutation Test Windows
all_perm_tests = fieldnames(perm);
for i = 1:length(all_perm_tests)
    test_name = all_perm_tests{i};
    tmp = find(perm.(test_name)(1, :) < p_val);
    id = tmp(consec_idx(tmp, thres));
    sig_windows.perm.(test_name) = find_sig_windows(id, time);
end

% Bootstrapping Test Windows
all_btsrp_tests = fieldnames(btsrp);
for i = 1:length(all_btsrp_tests)
    test_name = all_btsrp_tests{i};
    % Positive (signal > baseline)
    tmp_pos = find(btsrp.(test_name)(1, :) > 0);
    id_pos = tmp_pos(consec_idx(tmp_pos, thres));
    sig_windows.btsrp.(test_name).positive = find_sig_windows(id_pos, time);
    % Negative (signal < baseline)
    tmp_neg = find(btsrp.(test_name)(2, :) < 0);
    id_neg = tmp_neg(consec_idx(tmp_neg, thres));
    sig_windows.btsrp.(test_name).negative = find_sig_windows(id_neg, time);
end
fprintf('Significance window calculation complete.\n\n');
disp('SIGNIFICANCE WINDOWS (in seconds):');
disp(sig_windows);

% --- Plotting Section ---
fprintf('\nGenerating plots...\n');
% Plot 1: First Sessions - CSP vs CSM
if isfield(perm, 'first_csp_vs_csm')
    create_comparison_plot(CSP_first, CSM_first, btsrp.csp_first, btsrp.csm_first, perm.first_csp_vs_csm, ...
    {sprintf('CSP (n=%d)', size(CSP_first,1)), sprintf('CSM (n=%d)', size(CSM_first,1))}, ...
    {green, orange}, 'First Sessions: CSP vs CSM', time, p_val, thres);
end
% Plot 2: Final Sessions - CSP vs CSM
if isfield(perm, 'final_csp_vs_csm')
    create_comparison_plot(CSP_final, CSM_final, btsrp.csp_final, btsrp.csm_final, perm.final_csp_vs_csm, ...
    {sprintf('CSP (n=%d)', size(CSP_final,1)), sprintf('CSM (n=%d)', size(CSM_final,1))}, ...
    {green, orange}, 'Final Sessions: CSP vs CSM', time, p_val, thres);
end
% Plot 3: CSP - First vs Final Sessions
if isfield(perm, 'csp_first_vs_final')
    create_comparison_plot(CSP_first, CSP_final, btsrp.csp_first, btsrp.csp_final, perm.csp_first_vs_final, ...
    {sprintf('First Sessions (n=%d)', size(CSP_first,1)), sprintf('Final Sessions (n=%d)', size(CSP_final,1))}, ...
    {green, light_green}, 'CSP: First vs Final Sessions', time, p_val, thres);
end
% Plot 4: CSM - First vs Final Sessions
if isfield(perm, 'csm_first_vs_final')
    create_comparison_plot(CSM_first, CSM_final, btsrp.csm_first, btsrp.csm_final, perm.csm_first_vs_final, ...
    {sprintf('First Sessions (n=%d)', size(CSM_first,1)), sprintf('Final Sessions (n=%d)', size(CSM_final,1))}, ...
    {orange, light_orange}, 'CSM: First vs Final Sessions', time, p_val, thres);
end
fprintf('Plotting complete.\n');

%% ----------------- Plotting Function -----------------
function create_comparison_plot(data1, data2, btsrp1, btsrp2, perm_test_result, labels, colors, plot_title, time, p_val, thres)
    figure('Name', plot_title, 'NumberTitle', 'off');
    hold on;
    ylim([-2 3]); y_lim = ylim;
    
    datasets = {data1, data2};
    for i = 1:2
        mean_trace = mean(datasets{i}, 1);
        sem_trace = std(datasets{i}, 0, 1) / sqrt(size(datasets{i}, 1));
        plot(time, mean_trace, 'Color', colors{i}, 'LineWidth', 1.5);
        jbfill(time, mean_trace - sem_trace, mean_trace + sem_trace, colors{i}, 'none', 0, 0.2);
    end
    
    sig_y_perm = y_lim(1) + 0.5;
    tmp = find(perm_test_result(1, :) < p_val);
    id = tmp(consec_idx(tmp, thres));
    if ~isempty(id), plot(time(id), sig_y_perm * ones(size(id)), 's', 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'Color', 'k'); end
    
    all_btsrp_results = {btsrp1, btsrp2};
    boot_y_offsets = [y_lim(1) + 0.3, y_lim(1) + 0.1];
    for i = 1:2
        tmp_pos = find(all_btsrp_results{i}(1, :) > 0); id_pos = tmp_pos(consec_idx(tmp_pos, thres));
        tmp_neg = find(all_btsrp_results{i}(2, :) < 0); id_neg = tmp_neg(consec_idx(tmp_neg, thres));
        all_sig = [id_pos, id_neg];
        if ~isempty(all_sig), plot(time(all_sig), boot_y_offsets(i) * ones(size(all_sig)), 's', 'MarkerSize', 5, 'MarkerFaceColor', colors{i}, 'Color', colors{i}, 'HandleVisibility', 'off'); end
    end
    
    legend(labels, 'Location', 'northwest');
    title(plot_title); xlabel('Time from Cue Onset (s)'); ylabel('Z-Score');
    line([0 0], y_lim, 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
    line([10 10], y_lim, 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
    line(xlim, [0 0], 'Color', 'k', 'LineStyle', '-');
    box on; hold off;
end

%% ----------------- Helper Function for Significance Windows -----------------
function windows = find_sig_windows(indices, time_vector)
    % Takes a vector of significant indices and returns a struct array 
    % with the start and end times of each contiguous window.
    
    windows = struct('start_time_s', {}, 'end_time_s', {});
    if isempty(indices)
        return; % Return empty struct if no significant points
    end
    
    % Find breaks in consecutive indices
    breaks = find(diff(indices) > 1);
    
    % Define start and end points of each window
    start_indices = [indices(1), indices(breaks + 1)];
    end_indices = [indices(breaks), indices(end)];
    
    % Convert indices to times and store in the output struct
    for i = 1:length(start_indices)
        windows(i).start_time_s = time_vector(start_indices(i));
        windows(i).end_time_s = time_vector(end_indices(i));
    end
end