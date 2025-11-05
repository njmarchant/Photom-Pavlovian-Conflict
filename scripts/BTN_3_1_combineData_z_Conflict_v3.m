%% Nathan Marchant May 2025
% Written for Pavlovian conflict task
% --- MODIFIED: Now handles 'early' and 'late' session analysis ---
% --- MODIFIED: Plots early and late sessions separately ---
% --- MODIFIED: Plots a direct comparison of early vs. late for each cue type ---

%% --- Initialization and Setup ---
clear all;
close all;

% --- ‼️ ACTION REQUIRED: Define your sessions here ‼️ ---
% Add the session names (e.g., 'conflict_1') from 'sesdat.session' to the appropriate list.
early_sessions = {'early'}; %{'conflict_1', 'conflict_2'};%, 'conflict_3', 'conflict_4'};
late_sessions  = {'late'}; %{'conflict_11', 'conflict_12', 'conflict_13', 'conflict_14'};

% --- Define the folder containing your extracted .mat files ---
dataFolder = 'C:\Photometry\PavConf\Conflict_Extracted\Pavlovian conditioning conflict'; % <--- SET YOUR FOLDER PATH HERE
files = dir(fullfile(dataFolder, '*.mat'));
fprintf('Found %d data files in: %s\n', length(files), dataFolder);

% --- Initialize data matrices for Early and Late sessions ---
CSP_early = []; CSM_early = []; CSPun_early = [];
CSP_late  = []; CSM_late  = []; CSPun_late  = [];
CSP_early_source = {}; CSM_early_source = {}; CSPun_early_source = {};
CSP_late_source  = {}; CSM_late_source  = {}; CSPun_late_source  = {};

%% --- Data Extraction and Sorting Loop ---
for k = 1:length(files)
    currentFile = fullfile(dataFolder, files(k).name);
    load(currentFile, 'sesdat'); % Load only the 'sesdat' structure
    
    % Check if the loaded file contains the necessary data structure
    if exist('sesdat', 'var') && isfield(sesdat, 'session') && isfield(sesdat, 'traces_z') && ~isempty(sesdat.traces_z)
        if ~ismissing(sesdat.session_conf) 
            session_name = sesdat.session_conf;
            fprintf('Processing: %s (Session: %s)\n', files(k).name, session_name);
        else
           session_name = sesdat.session_conf;
            fprintf('Processing: %s (Session: %s)\n', files(k).name, 'empty');
        end
        
        % Find original row numbers for each trial type
        all_rows = (1:size(sesdat.traces_z, 1))';
        csp_rows = all_rows(sesdat.traces_z(:, 2) == 1);
        csm_rows = all_rows(sesdat.traces_z(:, 2) == 2);
        cspun_rows = all_rows(sesdat.traces_z(:, 2) == 3);
        
        % --- Sort data into EARLY or LATE based on session name ---
        if ismember(session_name, early_sessions)
            if ~isempty(csp_rows);   CSP_early = [CSP_early; sesdat.traces_z(csp_rows, 5:end)]; for r = 1:length(csp_rows); CSP_early_source{end+1,1} = sprintf('%s, row %d', files(k).name, csp_rows(r)); end; end
            if ~isempty(csm_rows);   CSM_early = [CSM_early; sesdat.traces_z(csm_rows, 5:end)]; for r = 1:length(csm_rows); CSM_early_source{end+1,1} = sprintf('%s, row %d', files(k).name, csm_rows(r)); end; end
            if ~isempty(cspun_rows); CSPun_early = [CSPun_early; sesdat.traces_z(cspun_rows, 5:end)]; for r = 1:length(cspun_rows); CSPun_early_source{end+1,1} = sprintf('%s, row %d', files(k).name, cspun_rows(r)); end; end
        elseif ismember(session_name, late_sessions)
            if ~isempty(csp_rows);   CSP_late = [CSP_late; sesdat.traces_z(csp_rows, 5:end)]; for r = 1:length(csp_rows); CSP_late_source{end+1,1} = sprintf('%s, row %d', files(k).name, csp_rows(r)); end; end
            if ~isempty(csm_rows);   CSM_late = [CSM_late; sesdat.traces_z(csm_rows, 5:end)]; for r = 1:length(csm_rows); CSM_late_source{end+1,1} = sprintf('%s, row %d', files(k).name, csm_rows(r)); end; end
            if ~isempty(cspun_rows); CSPun_late = [CSPun_late; sesdat.traces_z(cspun_rows, 5:end)]; for r = 1:length(cspun_rows); CSPun_late_source{end+1,1} = sprintf('%s, row %d', files(k).name, cspun_rows(r)); end; end
        else
            fprintf('  -> Skipping session "%s". Not in early/late lists.\n');
        end
    else
        fprintf('  -> Skipping file %s. "sesdat" structure invalid or empty.\n', files(k).name);
    end
    clear sesdat; % Clear sesdat to prevent data carry-over
end
fprintf('\nData collation complete.\n');

%% --- Data Cleaning ---

[CSP_early, CSP_early_source] = clean_data(CSP_early, CSP_early_source, 'CSP_early');
[CSM_early, CSM_early_source] = clean_data(CSM_early, CSM_early_source, 'CSM_early');
[CSPun_early, CSPun_early_source] = clean_data(CSPun_early, CSPun_early_source, 'CSPun_early');
[CSP_late, CSP_late_source] = clean_data(CSP_late, CSP_late_source, 'CSP_late');
[CSM_late, CSM_late_source] = clean_data(CSM_late, CSM_late_source, 'CSM_late');
[CSPun_late, CSPun_late_source] = clean_data(CSPun_late, CSPun_late_source, 'CSPun_late');

fprintf('\nData cleaning complete.\n');
fprintf('Total Early Trials (CSP/CSM/CSPun): %d / %d / %d\n', size(CSP_early, 1), size(CSM_early, 1), size(CSPun_early, 1));
fprintf('Total Late Trials (CSP/CSM/CSPun):  %d / %d / %d\n', size(CSP_late, 1), size(CSM_late, 1), size(CSPun_late, 1));

if isempty(CSP_early) && isempty(CSM_early) && isempty(CSPun_early) && isempty(CSP_late) && isempty(CSM_late) && isempty(CSPun_late)
    error('No data was loaded or all data was removed. Check file paths, session names, and data content.');
end

%% --- Plotting and Stats Configuration ---
% Colors (RGB values / 255)
green  = [0.10,0.60,0.00];
orange = [1.00,0.45,0.00];
red    = [1.0,0,0];
black1 = [0.75,0.75,0.75];
black2 = [0.5,0.5,0.5];
black3 = [0.25,0.25,0.25];


%% --- Plot 1: Early Sessions Summary ---
if ~isempty(CSP_early)
    time_axis = linspace(-10, 40, size(CSP_early, 2));
    datasets_early = {CSP_early, CSM_early, CSPun_early};
    labels_early = {sprintf('CSP (n=%d)', size(CSP_early, 1)), ...
                    sprintf('CSM (n=%d)', size(CSM_early, 1)), ...
                    sprintf('CSPun (n=%d)', size(CSPun_early, 1))};
    colors_early = {green, orange, red};
    plot_session_summary(time_axis, datasets_early, labels_early, colors_early, 'Early Sessions');
end

%% --- Plot 2: Late Sessions Summary ---
if ~isempty(CSP_late)
    time_axis = linspace(-10, 40, size(CSP_late, 2));
    datasets_late = {CSP_late, CSM_late, CSPun_late};
    labels_late = {sprintf('CSP (n=%d)', size(CSP_late, 1)), ...
                   sprintf('CSM (n=%d)', size(CSM_late, 1)), ...
                   sprintf('CSPun (n=%d)', size(CSPun_late, 1))};
    colors_late = {green, orange, red};
    plot_session_summary(time_axis, datasets_late, labels_late, colors_late, 'Late Sessions');
end

%% --- Plot 3: Early vs. Late Comparison for Each Cue Type ---
cue_types_early = {CSP_early, CSM_early, CSPun_early};
cue_types_late  = {CSP_late,  CSM_late,  CSPun_late};
cue_names       = {'CSP', 'CSM', 'CS Pun'};
cue_colors      = {green, orange, red};

for i = 1:length(cue_names)
    data_e = cue_types_early{i};
    data_l = cue_types_late{i};
    
    if isempty(data_e) || isempty(data_l)
        fprintf('Skipping %s Early vs. Late plot: Not enough data.\n', cue_names{i});
        continue;
    end
    
    figure('Name', sprintf('%s: Early vs. Late', cue_names{i}));
    hold on;
    
    time = linspace(-10, 40, size(data_e, 2));
    p = 0.01; thres = 8;
    
    % --- Stats: Permutation test between early and late ---
    [perm_e_l, ~] = permTest_array(data_e, data_l, 1000);
    
    % --- Plotting ---
    % Plot Late data (lighter shade)
    plot(time, mean(data_l, 1), 'Color', cue_colors{i}*0.6 + 0.4, 'LineWidth', 1.5);
    jbfill(time, (mean(data_l,1) - std(data_l,0,1)/sqrt(size(data_l,1))), ...
           (mean(data_l,1) + std(data_l,0,1)/sqrt(size(data_l,1))), ...
           cue_colors{i}*0.6 + 0.4, 'none', 0, 0.15);
           
    % Plot Early data (original color)
    plot(time, mean(data_e, 1), 'Color', cue_colors{i}, 'LineWidth', 1.5);
    jbfill(time, (mean(data_e,1) - std(data_e,0,1)/sqrt(size(data_e,1))), ...
           (mean(data_e,1) + std(data_e,0,1)/sqrt(size(data_e,1))), ...
           cue_colors{i}, 'none', 0, 0.15);
           
    % --- Legend and Labels ---
    legend({sprintf('%s Late (n=%d)', cue_names{i}, size(data_l, 1)), ...
            sprintf('%s Early (n=%d)', cue_names{i}, size(data_e, 1))}, ...
            'Location', 'northwest');
            
    % --- Plot Permutation Test Markers ---
    ax = gca; yLimits = ax.YLim; yMin = yLimits(1);
    tmp = find(perm_e_l(1,:) < p);
    id = tmp(consec_idx(tmp, thres));
    plot(time(id), (yMin - 0.2) * ones(size(id)), 's', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'Color', 'k');
    
    % --- Final Touches ---
    title(sprintf('%s: Early vs. Late', cue_names{i}));
    xlabel('Time (s)'); ylabel('Z-Scored Activity');
    line(ax.XLim, [0, 0], 'Color', 'k', 'LineStyle', '--');
    line([0, 0], ax.YLim, 'Color', 'k', 'LineStyle', '--');
    line([20, 20], ax.YLim, 'Color', 'k', 'LineStyle', '--');
    hold off;
end

%% --- Helper Functions ---
% ALL FUNCTION DEFINITIONS MUST BE AT THE END OF THE SCRIPT

function [data, source] = clean_data(data, source, name)
    % This helper function cleans a data matrix and its corresponding source list
    if isempty(data)
        return; % Do nothing if data is already empty
    end
    % Remove trials with extreme mean values
    mean_vals = mean(data, 2);
    outlier_idx = mean_vals > 20 | mean_vals < -20;
    if any(outlier_idx)
        fprintf('Found and removed %d outlier rows from %s.\n', sum(outlier_idx), name);
        data(outlier_idx, :) = [];
        if ~isempty(source); source(outlier_idx) = []; end
    end
    % Remove rows with any NaN values
    nan_idx = any(isnan(data), 2);
    if any(nan_idx)
        fprintf('Found and removed %d NaN rows from %s.\n', sum(nan_idx), name);
        data(nan_idx, :) = [];
        if ~isempty(source); source(nan_idx) = []; end
    end
end

function plot_session_summary(time, datasets, labels, colors, title_str)
    % This helper function plots a figure with 3 datasets (e.g., CSP, CSM, CSPun)
    % Colors (RGB values / 255)
    % green  = [0.10,0.60,0.00];
    % orange = [1.00,0.45,0.00];
    % red    = [1.0,0,0];
    black1 = [0.75,0.75,0.75];
    black2 = [0.5,0.5,0.5];
    black3 = [0.25,0.25,0.25];
    
    if all(cellfun('isempty', datasets))
        fprintf('Skipping "%s" plot: No data available.\n', title_str);
        return;
    end
    
    figure('Name', title_str);
    hold on;
    
    % --- Stats ---
    p = 0.01; thres = 8;
    btsrp = struct(); perm = struct();
    if ~isempty(datasets{1}); tmp = bootstrap_data(datasets{1}, 5000, 0.001); btsrp.d1 = CIadjust(tmp(1,:),tmp(2,:),tmp,size(datasets{1}, 1),2); end
    if ~isempty(datasets{2}); tmp = bootstrap_data(datasets{2}, 5000, 0.001); btsrp.d2 = CIadjust(tmp(1,:),tmp(2,:),tmp,size(datasets{2}, 1),2); end
    if ~isempty(datasets{3}); tmp = bootstrap_data(datasets{3}, 5000, 0.001); btsrp.d3 = CIadjust(tmp(1,:),tmp(2,:),tmp,size(datasets{3}, 1),2); end
    if ~isempty(datasets{1}) && ~isempty(datasets{2}); [perm.d1_d2, ~] = permTest_array(datasets{1}, datasets{2}, 1000); end
    if ~isempty(datasets{2}) && ~isempty(datasets{3}); [perm.d2_d3, ~] = permTest_array(datasets{2}, datasets{3}, 1000); end
    if ~isempty(datasets{1}) && ~isempty(datasets{3}); [perm.d1_d3, ~] = permTest_array(datasets{1}, datasets{3}, 1000); end
    
    % --- Plotting ---
    handles = gobjects(1, length(datasets));
    for i = 1:length(datasets)
        if ~isempty(datasets{i})
            handles(i) = plot(time, mean(datasets{i}, 1), 'Color', colors{i}, 'LineWidth', 1.5);
            jbfill(time, (mean(datasets{i}) - std(datasets{i},[],1)/sqrt(size(datasets{i},1))), ...
                   (mean(datasets{i}) + std(datasets{i},[],1)/sqrt(size(datasets{i},1))), colors{i}, 'none', 0, 0.2);
        end
    end
    legend(handles(isgraphics(handles)), labels(isgraphics(handles)), 'Location', 'northwest');
    
    % --- Stats Markers ---
    ax = gca; yLimits = ax.YLim; yMin = yLimits(1);
    offsetsperm = [-1.0, -1.2, -1.4]; markersperm = {black1, black2, black3};
    if isfield(perm, 'd1_d2'); tmp = find(perm.d1_d2(1,:) < p); id = tmp(consec_idx(tmp, thres)); plot(time(id), (yMin + offsetsperm(1)) * ones(size(id)), 's', 'MarkerSize', 7, 'MarkerFaceColor', markersperm{1}, 'Color', markersperm{1}); end
    if isfield(perm, 'd2_d3'); tmp = find(perm.d2_d3(1,:) < p); id = tmp(consec_idx(tmp, thres)); plot(time(id), (yMin + offsetsperm(2)) * ones(size(id)), 's', 'MarkerSize', 7, 'MarkerFaceColor', markersperm{2}, 'Color', markersperm{2}); end
    if isfield(perm, 'd1_d3'); tmp = find(perm.d1_d3(1,:) < p); id = tmp(consec_idx(tmp, thres)); plot(time(id), (yMin + offsetsperm(3)) * ones(size(id)), 's', 'MarkerSize', 7, 'MarkerFaceColor', markersperm{3}, 'Color', markersperm{3}); end
    
    % --- Final Touches ---
    line(ax.XLim, [0, 0], 'Color', 'k', 'LineStyle', '--');
    line([0, 0], ax.YLim, 'Color', 'k', 'LineStyle', '--');
    line([20, 20], ax.YLim, 'Color', 'k', 'LineStyle', '--');
    title(title_str); xlabel('Time (s)'); ylabel('Z-Scored Activity');
    hold off;
end