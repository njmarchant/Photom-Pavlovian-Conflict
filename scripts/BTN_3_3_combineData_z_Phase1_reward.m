%% Nathan Marchant May 2025
% Written for Pavlovian conflict task
% Phase 1 - reward
% Includes: automatic calculation of significance windows ---

% --- Initialization ---
clear all;
close all; % Good practice to close figures

dataFile = 'C:\Photometry\PavConf\IA_EtOH07AllData.mat'; % Make sure this is the correct file
load(dataFile); % This should load the 'alldata' structure

% Define the phase you want to analyze
PHASE_TO_ANALYZE = 'Reward'; % <--- SET YOUR DESIRED PHASE HERE ('Reward', 'Conflict', 'Discrimination')

% Define how CSP entries are grouped for early and late
early_csp_indices = [1, 2];
late_csp_indices = [3, 4];

% Initialize variables
fprintf('Loading data from: %s\n', dataFile);
if ~exist('alldata', 'var')
    error('"alldata" structure not found in the loaded file. Check file content or name.');
end
fprintf('Data loaded. Processing for phase: %s\n', PHASE_TO_ANALYZE);

EarlyCSP = []; 
LateCSP = [];
ALL_CSP = [];

for i = 1:length(alldata)
     currentRow = alldata(i);
     if isfield(currentRow, 'phase') && strcmp(currentRow.phase, PHASE_TO_ANALYZE)
         early_csp_data_1 = currentRow.csp{early_csp_indices(1)};
         early_csp_data_2 = currentRow.csp{early_csp_indices(2)};
                
         EarlyCSP = [EarlyCSP; early_csp_data_1; early_csp_data_2];
         
         late_csp_data_1 = currentRow.csp{late_csp_indices(1)};
         late_csp_data_2 = currentRow.csp{late_csp_indices(2)};
                
         LateCSP = [LateCSP; late_csp_data_1; late_csp_data_2];
     end
end 
ALL_CSP = [EarlyCSP; LateCSP]; % Combine after loops for clarity

%% 
%___________ colours used for the plots (RGB values / 255)______________
purp = [0.36,0.32,0.64]; blue = [0.00,0.26,1.00]; azure = [0.00,0.62,1.00];
red = [1.0,0,0]; pink = [1.00,0.50,0.50]; orange = [1.00,0.45,0.00];
light_orange = [1.00,0.87,0.60]; green = [0.10,0.60,0.00]; light_green = [0.00,0.50,0.50];
yellow = [1.00,0.88,0.00]; black1 = [0.75,0.75,0.75]; black2 = [0.5,0.5,0.5];
black3 = [0.25,0.25,0.25]; black4 = [0,0,0];

%___________ Dimensions and labels used for the plots ______________
time = linspace(-10, 40, size(EarlyCSP, 2));
Rew_Early_completed = num2str(size(EarlyCSP,1));
Rew_Early_Trial_label = ['CS+ Early (', Rew_Early_completed,')'];
Rew_Late_completed = num2str(size(LateCSP,1));
Rew_Late_Trial_label = ['CS+ Late (', Rew_Late_completed,')']; % Corrected label
Rew_All_completed = num2str(size(ALL_CSP,1));
Rew_All_Trial_label = ['CS+ All (', Rew_All_completed,')'];

%------STATS -----------------------------------------------------------------------
fprintf('\nPerforming statistical tests...\n');
% Bootstrap
btsrp = struct();
tmp = bootstrap_data(EarlyCSP, 5000, 0.001);
btsrp.rewE = CIadjust(tmp(1,:),tmp(2,:),tmp,size(EarlyCSP, 1),2);
tmp = bootstrap_data(LateCSP, 5000, 0.001);
btsrp.rewL = CIadjust(tmp(1,:),tmp(2,:),tmp,size(LateCSP, 1),2);
tmp = bootstrap_data(ALL_CSP, 5000, 0.001);
btsrp.allcsp = CIadjust(tmp(1,:),tmp(2,:),tmp,size(ALL_CSP, 1),2);

% Permutation tests
perm = struct();
[perm.rewE_rewL, ~] = permTest_array(EarlyCSP, LateCSP, 1000);
fprintf('Statistical tests complete.\n');

% --- Calculate and Store Significance Windows ---
fprintf('\nCalculating significance windows...\n');
% Define Significance Parameters
p = 0.01;    % p-value for permutation test
thres = 8;   % Consecutive data points threshold

sig_windows = struct();
% NOTE: Assumes 'consec_idx' is in your MATLAB path.

% Permutation Test Windows
all_perm_tests = fieldnames(perm);
for i = 1:length(all_perm_tests)
    test_name = all_perm_tests{i};
    tmp = find(perm.(test_name)(1, :) < p);
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


%-----------------------------------------------------Plot Reward (Early vs Late) -----------------------------------------------------------------------
fprintf('\nGenerating plots...\n');
figure('Name', 'Reward: Early vs Late');
datasets = {EarlyCSP, LateCSP};
labels = {Rew_Early_Trial_label, Rew_Late_Trial_label};
colors = {green, blue};
boottests = {'rewE', 'rewL'};
offsetsboot = [0, -0.1, -0.2, -0.3];
permtests = {'rewE_rewL'};
offsetsperm = -0.5;
markersperm = {yellow, black2, black3};

handles = zeros(1, numel(datasets));
for i = 1:length(colors)
    hold on
    handles(i) = plot(time, mean(datasets{i}, 1), 'Color', colors{i}, 'LineWidth', 1.5);
    jbfill(time, (mean(datasets{i}) - std(datasets{i}, 0, 1) / sqrt(size(datasets{i}, 1))), ...
        (std(datasets{i}, 0, 1) / sqrt(size(datasets{i}, 1)) + mean(datasets{i})), colors{i}, 'none', 0, 0.2);
end
legend(handles, labels, 'Location', 'northwest');

ax = gca;
yLimits = ax.YLim;
yMin = yLimits(1);

% Plotting markers for bootstrapping
for i = 1:length(boottests)
    % Positive
    tmp_pos = find(btsrp.(boottests{i})(1,:) > 0);
    id_pos = tmp_pos(consec_idx(tmp_pos, thres));
    plot(time(id_pos), (yMin + offsetsboot(i*2-1) * ones(size(time(id_pos)))), 's', 'MarkerSize', 7, 'MarkerFaceColor', colors{i}, 'Color', colors{i}, 'HandleVisibility', 'off');
    % Negative
    tmp_neg = find(btsrp.(boottests{i})(2,:) < 0);
    id_neg = tmp_neg(consec_idx(tmp_neg, thres));
    plot(time(id_neg), (yMin + offsetsboot(i*2) * ones(size(time(id_neg)))), 's', 'MarkerSize', 7, 'MarkerFaceColor', colors{i}, 'Color', colors{i}, 'HandleVisibility', 'off');
end

% Plotting permutation tests
for i = 1:length(permtests)
    tmp = find(perm.(permtests{i})(1, :) < p);
    id = tmp(consec_idx(tmp, thres));
    plot(time(id), (yMin + offsetsperm(i)) * ones(size(time(id))), 's', 'MarkerSize', 7, 'MarkerFaceColor', markersperm{i}, 'Color', markersperm{i}, 'HandleVisibility', 'off');
end

% Adding zero reference lines and title
yLimits = ax.YLim;
line(xlim, [0, 0], 'Color', 'black', 'HandleVisibility', 'off');
line([0, 0], yLimits, 'Color', 'black', 'HandleVisibility', 'off');
line([20, 20], yLimits, 'Color', 'black', 'HandleVisibility', 'off');
title('CS+ Phase 1: Early vs Late');  

%-----------------------------------------------------Plot Reward (All Trials) -----------------------------------------------------------------------
figure('Name', 'Reward: All Trials');
datasets = {ALL_CSP};
labels = {Rew_All_Trial_label};
colors = {green};
boottests = {'allcsp'};
offsetsboot = [0, -0.1];

hold on;
plot(time, mean(datasets{1}, 1), 'Color', colors{1}, 'LineWidth', 1.5);
jbfill(time, (mean(datasets{1}) - std(datasets{1}, 0, 1) / sqrt(size(datasets{1}, 1))), ...
    (std(datasets{1}, 0, 1) / sqrt(size(datasets{1}, 1)) + mean(datasets{1})), colors{1}, 'none', 0, 0.2);
legend(labels, 'Location', 'northwest');

ax = gca;
yLimits = ax.YLim;
yMin = yLimits(1);

% Plotting markers for bootstrapping
% Positive
tmp_pos = find(btsrp.allcsp(1,:) > 0);
id_pos = tmp_pos(consec_idx(tmp_pos, thres));
plot(time(id_pos), (yMin + offsetsboot(1) * ones(size(time(id_pos)))), 's', 'MarkerSize', 7, 'MarkerFaceColor', colors{1}, 'Color', colors{1}, 'HandleVisibility', 'off');
% Negative
tmp_neg = find(btsrp.allcsp(2,:) < 0);
id_neg = tmp_neg(consec_idx(tmp_neg, thres));
plot(time(id_neg), (yMin + offsetsboot(2) * ones(size(time(id_neg)))), 's', 'MarkerSize', 7, 'MarkerFaceColor', colors{1}, 'Color', colors{1}, 'HandleVisibility', 'off');

% Adding zero reference lines and title
yLimits = ax.YLim;
line(xlim, [0, 0], 'Color', 'black', 'HandleVisibility', 'off');
line([0, 0], yLimits, 'Color', 'black', 'HandleVisibility', 'off');
line([20, 20], yLimits, 'Color', 'black', 'HandleVisibility', 'off');
title('CS+ All Trials');
fprintf('Plotting complete.\n');

%% ----------------- Helper Function for Signficance Windows -----------------
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