%% Nathan Marchant May 2025
% Written for Pavlovian conflict task


% --- Initialization ---
clear all;
close all; % Good practice to close figures

% --- Define session groups ---
% List the session names to pool for the "early" and "late" groups.
early_sessions = {'shock_1', 'shock_2'};
late_sessions = {'shock_3', 'shock_4'};

% --- Define the folder containing your extracted .mat files ---
dataFolder = 'C:\Photometry\PavConf\DrPhotom_Extracted\Alcohol conditioning V4 photometry_SHOCK'; % <--- SET YOUR FOLDER PATH HERE
files = dir(fullfile(dataFolder, '*.mat'));
fprintf('Found %d data files in: %s\n', length(files), dataFolder);

% --- Initialize data matrices for each condition and session group ---
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

%% ----------------- Plot 1: 'Early' sessions comparison -----------------
if ~isempty(CSP_no_shock_early) || ~isempty(CSP_shock_early) || ~isempty(CSM_early)
    figure('Name', 'Early Sessions', 'NumberTitle', 'off');
    hold on;
    
    % --- Stats ---
    if ~isempty(CSP_no_shock_early) && ~isempty(CSP_shock_early), [perm.early_1, ~] = permTest_array(CSP_no_shock_early, CSP_shock_early, 1000); else perm.early_1 = []; end
    if ~isempty(CSP_no_shock_early) && ~isempty(CSM_early), [perm.early_2, ~] = permTest_array(CSP_no_shock_early, CSM_early, 1000); else perm.early_2 = []; end
    if ~isempty(CSP_shock_early) && ~isempty(CSM_early), [perm.early_3, ~] = permTest_array(CSP_shock_early, CSM_early, 1000); else perm.early_3 = []; end

    % --- Plotting ---
    all_datasets = {CSP_no_shock_early, CSP_shock_early, CSM_early};
    all_labels = {sprintf('CS+ No Shock (n=%d)', size(CSP_no_shock_early,1)), sprintf('CS+ Shock (n=%d)', size(CSP_shock_early,1)), sprintf('CS- (n=%d)', size(CSM_early,1))};
    all_colors = {green, red, orange};
    
    plot_handles = [];
    legend_labels = {};
    for i = 1:length(all_datasets)
        if ~isempty(all_datasets{i})
            h = plot(time, mean(all_datasets{i}, 1), 'Color', all_colors{i}, 'LineWidth', 1.5);
            jbfill(time, (mean(all_datasets{i},1) - std(all_datasets{i},0,1)/sqrt(size(all_datasets{i},1))), (mean(all_datasets{i},1) + std(all_datasets{i},0,1)/sqrt(size(all_datasets{i},1))), all_colors{i}, 'none', 0.2);
            plot_handles = [plot_handles, h];
            legend_labels{end+1} = all_labels{i};
        end
    end
    
    if ~isempty(plot_handles), legend(plot_handles, legend_labels, 'Location', 'northwest'); end
    ax = gca; yLimits = ax.YLim; yMin = yLimits(1);
    
    % --- Significance Bars ---
    offsetsperm = [yMin + 0.05 * (yLimits(2)-yLimits(1)), yMin + 0.10 * (yLimits(2)-yLimits(1)), yMin + 0.15 * (yLimits(2)-yLimits(1))];
    permtests = {perm.early_1, perm.early_2, perm.early_3};
    markersperm = {black1, black2, black3};
    for i = 1:length(permtests)
        if ~isempty(permtests{i})
            tmp = find(permtests{i}(1, :) < p);
            id = tmp(consec_idx(tmp, thres));
            plot(time(id), offsetsperm(i) * ones(size(id)), 's', 'MarkerSize', 7, 'MarkerFaceColor', markersperm{i}, 'Color', markersperm{i}, 'HandleVisibility', 'off');
        end
    end
    
    % --- Final Touches ---
    title('Early Sessions'); xlabel('Time from Cue Onset (s)'); ylabel('Z-Score');
    line([0 0], ylim, 'Color', [0.5 0.5 0.5], 'LineStyle', '--'); line([10 10], ylim, 'Color', [0.5 0.5 0.5], 'LineStyle', '--'); line(xlim, [0 0], 'Color', 'k');
    box on; hold off;
end

%% ----------------- Plot 2: 'Late' sessions comparison ------------------
if ~isempty(CSP_no_shock_late) || ~isempty(CSP_shock_late) || ~isempty(CSM_late)
    figure('Name', 'Late Sessions', 'NumberTitle', 'off');
    hold on;
    
    % --- Stats ---
    if ~isempty(CSP_no_shock_late) && ~isempty(CSP_shock_late), [perm.late_1, ~] = permTest_array(CSP_no_shock_late, CSP_shock_late, 1000); else perm.late_1 = []; end
    if ~isempty(CSP_no_shock_late) && ~isempty(CSM_late), [perm.late_2, ~] = permTest_array(CSP_no_shock_late, CSM_late, 1000); else perm.late_2 = []; end
    if ~isempty(CSP_shock_late) && ~isempty(CSM_late), [perm.late_3, ~] = permTest_array(CSP_shock_late, CSM_late, 1000); else perm.late_3 = []; end

    % --- Plotting ---
    all_datasets = {CSP_no_shock_late, CSP_shock_late, CSM_late};
    all_labels = {sprintf('CS+ No Shock (n=%d)', size(CSP_no_shock_late,1)), sprintf('CS+ Shock (n=%d)', size(CSP_shock_late,1)), sprintf('CS- (n=%d)', size(CSM_late,1))};
    all_colors = {green, red, orange};
    
    plot_handles = [];
    legend_labels = {};
    for i = 1:length(all_datasets)
        if ~isempty(all_datasets{i})
            h = plot(time, mean(all_datasets{i}, 1), 'Color', all_colors{i}, 'LineWidth', 1.5);
            jbfill(time, (mean(all_datasets{i},1) - std(all_datasets{i},0,1)/sqrt(size(all_datasets{i},1))), (mean(all_datasets{i},1) + std(all_datasets{i},0,1)/sqrt(size(all_datasets{i},1))), all_colors{i}, 'none', 0.2);
            plot_handles = [plot_handles, h];
            legend_labels{end+1} = all_labels{i};
        end
    end
    
    if ~isempty(plot_handles), legend(plot_handles, legend_labels, 'Location', 'northwest'); end
    ax = gca; yLimits = ax.YLim; yMin = yLimits(1);
    
    % --- Significance Bars ---
    offsetsperm = [yMin + 0.05 * (yLimits(2)-yLimits(1)), yMin + 0.10 * (yLimits(2)-yLimits(1)), yMin + 0.15 * (yLimits(2)-yLimits(1))];
    permtests = {perm.late_1, perm.late_2, perm.late_3};
    markersperm = {black1, black2, black3};
    for i = 1:length(permtests)
        if ~isempty(permtests{i})
            tmp = find(permtests{i}(1, :) < p);
            id = tmp(consec_idx(tmp, thres));
            plot(time(id), offsetsperm(i) * ones(size(id)), 's', 'MarkerSize', 7, 'MarkerFaceColor', markersperm{i}, 'Color', markersperm{i}, 'HandleVisibility', 'off');
        end
    end
    
    % --- Final Touches ---
    title('Late Sessions'); xlabel('Time from Cue Onset (s)'); ylabel('Z-Score');
    line([0 0], ylim, 'Color', [0.5 0.5 0.5], 'LineStyle', '--'); line([10 10], ylim, 'Color', [0.5 0.5 0.5], 'LineStyle', '--'); line(xlim, [0 0], 'Color', 'k');
    box on; hold off;
end

%% ----------------- Plot 3: CS+ No Shock (Early vs Late) ----------------
if ~isempty(CSP_no_shock_early) || ~isempty(CSP_no_shock_late)
    figure('Name', 'CS+ No Shock: Early vs Late', 'NumberTitle', 'off');
    hold on;
    
    % --- Stats ---
    if ~isempty(CSP_no_shock_early) && ~isempty(CSP_no_shock_late), [perm.cspns_e_l, ~] = permTest_array(CSP_no_shock_early, CSP_no_shock_late, 1000); else perm.cspns_e_l = []; end

    % --- Plotting ---
    all_datasets = {CSP_no_shock_early, CSP_no_shock_late};
    all_labels = {sprintf('Early (n=%d)', size(CSP_no_shock_early,1)), sprintf('Late (n=%d)', size(CSP_no_shock_late,1))};
    all_colors = {green, light_green};
    
    plot_handles = [];
    legend_labels = {};
    for i = 1:length(all_datasets)
        if ~isempty(all_datasets{i})
            h = plot(time, mean(all_datasets{i}, 1), 'Color', all_colors{i}, 'LineWidth', 1.5);
            jbfill(time, (mean(all_datasets{i},1) - std(all_datasets{i},0,1)/sqrt(size(all_datasets{i},1))), (mean(all_datasets{i},1) + std(all_datasets{i},0,1)/sqrt(size(all_datasets{i},1))), all_colors{i}, 'none', 0.2);
            plot_handles = [plot_handles, h];
            legend_labels{end+1} = all_labels{i};
        end
    end
    
    if ~isempty(plot_handles), legend(plot_handles, legend_labels, 'Location', 'northwest'); end
    ax = gca; yLimits = ax.YLim; yMin = yLimits(1);
    
    % --- Significance Bar ---
    if ~isempty(perm.cspns_e_l)
        tmp = find(perm.cspns_e_l(1, :) < p);
        id = tmp(consec_idx(tmp, thres));
        plot(time(id), (yMin + 0.05 * (yLimits(2)-yLimits(1))) * ones(size(id)), 's', 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'Color', 'k', 'HandleVisibility', 'off');
    end
    
    % --- Final Touches ---
    title('CS+ No Shock: Early vs Late'); xlabel('Time from Cue Onset (s)'); ylabel('Z-Score');
    line([0 0], ylim, 'Color', [0.5 0.5 0.5], 'LineStyle', '--'); line([10 10], ylim, 'Color', [0.5 0.5 0.5], 'LineStyle', '--'); line(xlim, [0 0], 'Color', 'k');
    box on; hold off;
end

%% ----------------- Plot 4: CS+ Shock (Early vs Late) -------------------
if ~isempty(CSP_shock_early) || ~isempty(CSP_shock_late)
    figure('Name', 'CS+ Shock: Early vs Late', 'NumberTitle', 'off');
    hold on;
    
    % --- Stats ---
    if ~isempty(CSP_shock_early) && ~isempty(CSP_shock_late), [perm.csps_e_l, ~] = permTest_array(CSP_shock_early, CSP_shock_late, 1000); else perm.csps_e_l = []; end

    % --- Plotting ---
    all_datasets = {CSP_shock_early, CSP_shock_late};
    all_labels = {sprintf('Early (n=%d)', size(CSP_shock_early,1)), sprintf('Late (n=%d)', size(CSP_shock_late,1))};
    all_colors = {red, light_red};
    
    plot_handles = [];
    legend_labels = {};
    for i = 1:length(all_datasets)
        if ~isempty(all_datasets{i})
            h = plot(time, mean(all_datasets{i}, 1), 'Color', all_colors{i}, 'LineWidth', 1.5);
            jbfill(time, (mean(all_datasets{i},1) - std(all_datasets{i},0,1)/sqrt(size(all_datasets{i},1))), (mean(all_datasets{i},1) + std(all_datasets{i},0,1)/sqrt(size(all_datasets{i},1))), all_colors{i}, 'none', 0.2);
            plot_handles = [plot_handles, h];
            legend_labels{end+1} = all_labels{i};
        end
    end
    
    if ~isempty(plot_handles), legend(plot_handles, legend_labels, 'Location', 'northwest'); end
    ax = gca; yLimits = ax.YLim; yMin = yLimits(1);
    
    % --- Significance Bar ---
    if ~isempty(perm.csps_e_l)
        tmp = find(perm.csps_e_l(1, :) < p);
        id = tmp(consec_idx(tmp, thres));
        plot(time(id), (yMin + 0.05 * (yLimits(2)-yLimits(1))) * ones(size(id)), 's', 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'Color', 'k', 'HandleVisibility', 'off');
    end
    
    % --- Final Touches ---
    title('CS+ Shock: Early vs Late'); xlabel('Time from Cue Onset (s)'); ylabel('Z-Score');
    line([0 0], ylim, 'Color', [0.5 0.5 0.5], 'LineStyle', '--'); line([10 10], ylim, 'Color', [0.5 0.5 0.5], 'LineStyle', '--'); line(xlim, [0 0], 'Color', 'k');
    box on; hold off;
end

%% ----------------- Plot 5: CS- (Early vs Late) -------------------------
if ~isempty(CSM_early) || ~isempty(CSM_late)
    figure('Name', 'CS-: Early vs Late', 'NumberTitle', 'off');
    hold on;
    
    % --- Stats ---
    if ~isempty(CSM_early) && ~isempty(CSM_late), [perm.csm_e_l, ~] = permTest_array(CSM_early, CSM_late, 1000); else perm.csm_e_l = []; end

    % --- Plotting ---
    all_datasets = {CSM_early, CSM_late};
    all_labels = {sprintf('Early (n=%d)', size(CSM_early,1)), sprintf('Late (n=%d)', size(CSM_late,1))};
    all_colors = {orange, light_orange};
    
    plot_handles = [];
    legend_labels = {};
    for i = 1:length(all_datasets)
        if ~isempty(all_datasets{i})
            h = plot(time, mean(all_datasets{i}, 1), 'Color', all_colors{i}, 'LineWidth', 1.5);
            jbfill(time, (mean(all_datasets{i},1) - std(all_datasets{i},0,1)/sqrt(size(all_datasets{i},1))), (mean(all_datasets{i},1) + std(all_datasets{i},0,1)/sqrt(size(all_datasets{i},1))), all_colors{i}, 'none', 0.2);
            plot_handles = [plot_handles, h];
            legend_labels{end+1} = all_labels{i};
        end
    end
    
    if ~isempty(plot_handles), legend(plot_handles, legend_labels, 'Location', 'northwest'); end
    ax = gca; yLimits = ax.YLim; yMin = yLimits(1);
    
    % --- Significance Bar ---
    if ~isempty(perm.csm_e_l)
        tmp = find(perm.csm_e_l(1, :) < p);
        id = tmp(consec_idx(tmp, thres));
        plot(time(id), (yMin + 0.05 * (yLimits(2)-yLimits(1))) * ones(size(id)), 's', 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'Color', 'k', 'HandleVisibility', 'off');
    end
    
    % --- Final Touches ---
    title('CS-: Early vs Late'); xlabel('Time from Cue Onset (s)'); ylabel('Z-Score');
    line([0 0], ylim, 'Color', [0.5 0.5 0.5], 'LineStyle', '--'); line([10 10], ylim, 'Color', [0.5 0.5 0.5], 'LineStyle', '--'); line(xlim, [0 0], 'Color', 'k');
    box on; hold off;
end
