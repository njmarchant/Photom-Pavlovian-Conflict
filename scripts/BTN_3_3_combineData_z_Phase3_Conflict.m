%% Nathan Marchant May 2025
% Written for Pavlovian conflict task
% Phase 1 - reward

% --- Initialization ---
clear all;
close all; % Good practice to close figures
dataFile = 'C:\Photometry\PavConf\IA_EtOH07AllData.mat'; % Make sure this is the correct file
load(dataFile); % This should load the 'alldata' structure

% Define the phase you want to analyze
PHASE_TO_ANALYZE = 'Conflict'; % <--- SET YOUR DESIRED PHASE HERE ('Reward', 'Conflict', 'Discrimination')

% Define how CSP entries are grouped for early and late
% These are based on your previous clarification for 'Reward' phase
early_csp_indices = [1, 2];
late_csp_indices = [3, 4];
expected_csp_count_for_phase = 4; % For 'Reward' phase

% Initialize variables to store aggregated early and late trial data
EarlyCSP = [];
LateCSP = [];

fprintf('Loading data from: %s\n', dataFile);
if ~exist('alldata', 'var')
    error('"alldata" structure not found in the loaded file. Check file content or name.');
end
fprintf('Data loaded. Processing for phase: %s\n', PHASE_TO_ANALYZE);

% --- Data Extraction and Processing Loop ---
EarlyCSP = []; % Initialize outside
LateCSP = [];  % Initialize outside
EarlyCSM = []; % Initialize outside
LateCSM = [];  % Initialize outside
EarlyCSPun = []; % Initialize outside
LateCSPun = [];  % Initialize outside
for i = 1:length(alldata)
     currentRow = alldata(i);
     if isfield(currentRow, 'phase') && strcmp(currentRow.phase, PHASE_TO_ANALYZE)
         early_csp_data_1 = currentRow.csp{early_csp_indices(1)};
         early_csp_data_2 = currentRow.csp{early_csp_indices(2)};
         early_csm_data_1 = currentRow.csm{early_csp_indices(1)};
         early_csm_data_2 = currentRow.csm{early_csp_indices(2)};
         early_cspun_data_1 = currentRow.cspun{early_csp_indices(1)};
         early_cspun_data_2 = currentRow.cspun{early_csp_indices(2)};
                
         EarlyCSP = [EarlyCSP;early_csp_data_1; early_csp_data_2];
         EarlyCSM = [EarlyCSM;early_csm_data_1; early_csm_data_2];
         EarlyCSPun = [EarlyCSPun;early_cspun_data_1; early_cspun_data_2];

         late_csp_data_1 = currentRow.csp{late_csp_indices(1)};
         late_csp_data_2 = currentRow.csp{late_csp_indices(2)};
         late_csm_data_1 = currentRow.csm{late_csp_indices(1)};
         late_csm_data_2 = currentRow.csm{late_csp_indices(2)};
         late_cspun_data_1 = currentRow.cspun{late_csp_indices(1)};
         late_cspun_data_2 = currentRow.cspun{late_csp_indices(2)};
                
         LateCSP = [LateCSP;late_csp_data_1; late_csp_data_2];
         LateCSM = [LateCSM;late_csm_data_1; late_csm_data_2];
         LateCSPun = [LateCSPun;late_cspun_data_1; late_cspun_data_2];
            
     end
end


% --- Safety check after loop ---
if isempty(EarlyCSP) || isempty(LateCSP)
    warning('No data collected for phase: %s. Check PHASE_TO_ANALYZE and data structure.', PHASE_TO_ANALYZE);
    % You might want to return or error here if no data is found,
    % otherwise the plotting code will error.
    return; 
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
 
    time = linspace(-10, 40, size(EarlyCSP, 2));
    Rew_Early_completed = num2str(size(EarlyCSP,1));
    Rew_Early_Trial_label = ['CS+ Early (', Rew_Early_completed,')'];
    Rew_Late_completed = num2str(size(LateCSP,1));
    Rew_Late_Trial_label = ['CS+ Late (', Rew_Late_completed,')'];

    CSM_Early_completed = num2str(size(EarlyCSM,1));
    CSM_Early_Trial_label = ['CSM Early (', CSM_Early_completed,')'];
    CSM_Late_completed = num2str(size(LateCSM,1));
    CSM_Late_Trial_label = ['CS+ Late (', CSM_Late_completed,')'];

    CSPun_Early_completed = num2str(size(EarlyCSPun,1));
    CSPun_Early_Trial_label = ['CS+ Early (', CSPun_Early_completed,')'];
    CSPun_Late_completed = num2str(size(LateCSPun,1));
    CSPun_Late_Trial_label = ['CS+ Late (', CSPun_Late_completed,')'];
    

            %------STATS -----------------------------------------------------------------------

% %bootstrap 
tmp = bootstrap_data(EarlyCSP, 5000, 0.001);
btsrp.rewE = CIadjust(tmp(1,:),tmp(2,:),tmp,size(EarlyCSP, 1),2);
tmp = bootstrap_data(LateCSP, 5000, 0.001);
btsrp.rewL = CIadjust(tmp(1,:),tmp(2,:),tmp,size(LateCSP, 1),2);

tmp = bootstrap_data(EarlyCSM, 5000, 0.001);
btsrp.csmE = CIadjust(tmp(1,:),tmp(2,:),tmp,size(EarlyCSM, 1),2);
tmp = bootstrap_data(LateCSM, 5000, 0.001);
btsrp.csmL = CIadjust(tmp(1,:),tmp(2,:),tmp,size(LateCSM, 1),2);

tmp = bootstrap_data(EarlyCSPun, 5000, 0.001);
btsrp.cspunE = CIadjust(tmp(1,:),tmp(2,:),tmp,size(EarlyCSPun, 1),2);
tmp = bootstrap_data(LateCSPun, 5000, 0.001);
btsrp.cspunL = CIadjust(tmp(1,:),tmp(2,:),tmp,size(LateCSPun, 1),2);


% %permutation tests
[perm.rewE_rewL, ~] = permTest_array(EarlyCSP, LateCSP, 1000);
[perm.csmE_csmL, ~] = permTest_array(EarlyCSM, LateCSM, 1000);
[perm.cspunE_cspunL, ~] = permTest_array(EarlyCSPun, LateCSPun, 1000);


%-----------------------------------------------------Plot CSP  -----------------------------------------------------------------------
 
%Statistical parameters
p = 0.01;
thres = 8;

figure
datasets = {EarlyCSP, LateCSP};
labels = {Rew_Early_Trial_label, Rew_Late_Trial_label};
colors = {red, blue};

boottests = {'rewL', 'rewL'};
offsetsboot = [0, -0.1, -0.2, -0.3];
permtests = {'rewE_rewL'};
offsetsperm = -0.5;
 markersperm = {yellow, black2, black3};

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
line([20, 20], [yLimits(1), yLimits(2)], 'Color', 'black', 'HandleVisibility', 'off');
title('LateCSPun v LateCSP Phase 3');  


%-----------------------------------------------------Plot CS M-----------------------------------------------------------------------
 
%Statistical parameters
p = 0.01;
thres = 8;

figure
datasets = {EarlyCSM, LateCSM};
labels = {CSM_Early_Trial_label, CSM_Late_Trial_label};
colors = {green, azure};

boottests = {'csmE', 'csmL'};
offsetsboot = [0, -0.1, -0.2, -0.3];
permtests = {'csmE_csmL'};
offsetsperm = -0.5;
 markersperm = {yellow, black2, black3};

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
line([20, 20], [yLimits(1), yLimits(2)], 'Color', 'black', 'HandleVisibility', 'off');
title('CS M Phase 3'); 

%-----------------------------------------------------Plot CSPUN  -----------------------------------------------------------------------
 
%Statistical parameters
p = 0.01;
thres = 8;

figure
datasets = {EarlyCSPun, LateCSPun};
labels = {CSPun_Early_Trial_label, CSPun_Late_Trial_label};
colors = {red, orange};

boottests = {'cspunE', 'cspunL'};
offsetsboot = [0, -0.1, -0.2, -0.3];
permtests = {'cspunE_cspunL'};
offsetsperm = -0.5;
 markersperm = {yellow, black2, black3};

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
line([20, 20], [yLimits(1), yLimits(2)], 'Color', 'black', 'HandleVisibility', 'off');
title('CS PUN Phase 3'); 