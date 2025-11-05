%% Nathan Marchant May 2025
% Written for Pavlovian conflict task
% --- MODIFIED: Now loads and collates data from all .mat files in a folder ---
% --- MODIFIED: Added handling for NaN values ---
% --- MODIFIED: Added data source tracking and outlier removal based on mean value ---

% --- Initialization ---
clear all;
close all; % Good practice to close figures

% --- MODIFIED: Define the folder containing your extracted .mat files ---
dataFolder = 'C:\Photometry\PavConf\DrPhotom_Extracted\Alcohol conditioning V4 photometry'; % <--- SET YOUR FOLDER PATH HERE
files = dir(fullfile(dataFolder, '*.mat')); % Get a list of all .mat files
fprintf('Found %d data files in: %s\n', length(files), dataFolder);

% --- Data Extraction and Processing Loop ---
CSP = []; CSM = []; % Initialize data matrices
CSP_source = {}; CSM_source = {}; % Initialize cell arrays to track data origin

for k = 1:length(files)
    currentFile = fullfile(dataFolder, files(k).name);
    fprintf('Loading and processing: %s\n', files(k).name);
    
    load(currentFile); % This loads the 'sesdat' structure
    
    % Check if the loaded file contains the necessary data structure
    if exist('sesdat', 'var') && isfield(sesdat, 'traces_z') && ~isempty(sesdat.traces_z)
        
        % Find original row numbers for each trial type
        all_rows = (1:size(sesdat.traces_z, 1))';
        csp_orig_rows = all_rows(sesdat.traces_z(:, 2) == 1);
        csm_orig_rows = all_rows(sesdat.traces_z(:, 2) == 2);
                
        % Append the trace data and its source information
        if ~isempty(csp_orig_rows)
            CSP = [CSP; sesdat.traces_z(csp_orig_rows, 5:end)];
            for r = 1:length(csp_orig_rows)
                CSP_source{end+1, 1} = sprintf('%s, original row %d', files(k).name, csp_orig_rows(r));
            end
        end
        if ~isempty(csm_orig_rows)
            CSM = [CSM; sesdat.traces_z(csm_orig_rows, 5:end)];
            for r = 1:length(csm_orig_rows)
                CSM_source{end+1, 1} = sprintf('%s, original row %d', files(k).name, csm_orig_rows(r));
            end
        end
        
    else
        fprintf('  -> Skipping file. "sesdat.traces_z" not found or is empty.\n');
    end
    
    clear sesdat; % Clear sesdat to prevent data carry-over
end

fprintf('\nData collation complete.\n');


% --- MODIFIED: Remove trials with extreme mean values ---
mean_csp = mean(CSP, 2);
outlier_idx_csp = mean_csp > 20 | mean_csp < -20;
if any(outlier_idx_csp)
    fprintf('Found and removed %d rows from CSP with mean outside +/- 20.\n', sum(outlier_idx_csp));
    CSP(outlier_idx_csp, :) = [];
    CSP_source(outlier_idx_csp) = []; % Keep source info synchronized
end

mean_csm = mean(CSM, 2);
outlier_idx_csm = mean_csm > 20 | mean_csm < -20;
if any(outlier_idx_csm)
    fprintf('Found and removed %d rows from CSM with mean outside +/- 20.\n', sum(outlier_idx_csm));
    CSM(outlier_idx_csm, :) = [];
    CSM_source(outlier_idx_csm) = [];
end



% --- MODIFIED: Handle NaN values ---
nan_rows_csp = any(isnan(CSP), 2);
if any(nan_rows_csp)
    fprintf('Found and removed %d rows with NaN values from CSP.\n', sum(nan_rows_csp));
    CSP(nan_rows_csp, :) = [];
end

nan_rows_csm = any(isnan(CSM), 2);
if any(nan_rows_csm)
    fprintf('Found and removed %d rows with NaN values from CSM.\n', sum(nan_rows_csm));
    CSM(nan_rows_csm, :) = [];
end




fprintf('\nData cleaning complete.\n');
fprintf('Total CSP trials after cleaning: %d\n', size(CSP, 1));
fprintf('Total CSM trials after cleaning: %d\n', size(CSM, 1));


if isempty(CSP) && isempty(CSM) 
    error('No data was loaded or all data was removed. Check the folder path and the content of your .mat files.');
end
 %% 
 %___________ colours used for the plots (RGB values / 255)______________
    purp = [0.36,0.32,0.64];
    blue = [0.00,0.26,1.00];
    azure = [0.00,0.62,1.00];
    red = [1.0,0,0];
    pink = [1.00,0.50,0.50];
    orange = [1.00,0.45,0.00];
    light_orange = [1.00,0.87,0.60];
    green = [0.10,0.60,0.00];
    light_green = [0.00,0.50,0.50];
    yellow = [1.00,0.88,0.00];
    black1 = [0.75,0.75,0.75];
    black2 = [0.5,0.5,0.5];
    black3 = [0.25,0.25,0.25];
    black4 = [0,0,0];
%___________ Dimensions and labels used for the plots ______________
 
    time = linspace(-10, 20, size(CSP, 2));
    Rew_Early_completed = num2str(size(CSP,1));
    Rew_Early_Trial_label = ['CSP (', Rew_Early_completed,')'];
    CSM_Early_completed = num2str(size(CSM,1));
    CSM_Early_Trial_label = ['CSM (', CSM_Early_completed,')'];
  
    
    
    
%% ------STATS -----------------------------------------------------------------------
% %bootstrap 
tmp = bootstrap_data(CSP, 5000, 0.001);
btsrp.csp = CIadjust(tmp(1,:),tmp(2,:),tmp,size(CSP, 1),2);
tmp = bootstrap_data(CSM, 5000, 0.001);
btsrp.csm = CIadjust(tmp(1,:),tmp(2,:),tmp,size(CSM, 1),2);

% %permutation tests
[perm.csp_csm, ~] = permTest_array(CSP, CSM, 1000);

%% -----------------------------------------------------Plot CSP  -----------------------------------------------------------------------
 
%Statistical parameters
p = 0.01;
thres = 8;
figure
datasets = {CSP, CSM};
labels = {Rew_Early_Trial_label, CSM_Early_Trial_label};
colors = {green, orange};
boottests = {'csp', 'csm'};
offsetsboot = [0, -0.1, -0.2, -0.3, -0.4, -0.5];
permtests = {'csp_csm'};
offsetsperm = [-1.0, -1.2, -1.4, -1.6];
markersperm = {black1, black2, black3, black4};
handles = zeros(1, numel(datasets)*2); % Preallocate handles
for i = 1:length(colors)
    hold on
    handles(i) = plot(time, mean(datasets{i}, 1), 'Color', colors{i});
    jbfill(time, (mean(datasets{i}) - std(datasets{i}, 0, 1) / sqrt(size(datasets{i}, 1))), ...
        (std(datasets{i}, 0, 1) / sqrt(size(datasets{i}, 1)) + mean(datasets{i})), colors{i}, 'none', 0, 0.2);
    handles(i+3) = plot(NaN, NaN, 'Color', colors{i}); % Dummy handle for legend
end
legend(handles(1:length(labels)), labels, 'Location', 'northwest'); % Provide only the handles for mean lines
ax = gca;
yLimits = ax.YLim;
yMin = yLimits(1);
% Plotting markers for bootstrapping
for i = 1:2:length(boottests)*2
    tmp = find(btsrp.(boottests{(i+1)/2})(1,:) > 0);   % Find indices for values > 0 (i.e. signal higher than baseline)
    id = tmp(consec_idx(tmp, thres));
    plot(time(id), (yMin + offsetsboot(i) * ones(size(time(id), 2), 2)), 's', 'MarkerSize', 7, 'MarkerFaceColor', colors{(i+1)/2}, 'Color', colors{(i+1)/2}, 'HandleVisibility', 'off');
    clear tmp id
    tmp = find(btsrp.(boottests{(i+1)/2})(2,:) < 0); % Find indices for values < 0 (i.e. signal lower than baseline)
    id = tmp(consec_idx(tmp, thres));
    plot(time(id), (yMin + offsetsboot(i+1) * ones(size(time(id), 2), 2)), 's', 'MarkerSize', 7, 'MarkerFaceColor', colors{(i+1)/2}, 'Color', colors{(i+1)/2}, 'HandleVisibility', 'off');
end
% Plotting permutation tests
for i = 1:length(permtests)
    tmp = find(perm.(permtests{i})(1, :) < p);
    id = tmp(consec_idx(tmp, thres));
    plot(time(id), (yMin + offsetsperm(i)) * ones(size(time(id), 2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor', markersperm{i}, 'Color', markersperm{i}, 'HandleVisibility', 'off');
end
% Adding zero reference lines and title
yLimits = ax.YLim;
line([ax.XLim(1), ax.XLim(2)], [0, 0], 'Color', 'black', 'HandleVisibility', 'off');
line([0, 0], [yLimits(1), yLimits(2)], 'Color', 'black', 'HandleVisibility', 'off');
line([10, 10], [yLimits(1), yLimits(2)], 'Color', 'black', 'HandleVisibility', 'off');
title(' Phase 3');
