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




fprintf('Loading data from: %s\n', dataFile);
if ~exist('alldata', 'var')
    error('"alldata" structure not found in the loaded file. Check file content or name.');
end
fprintf('Data loaded. Processing for phase: %s\n', PHASE_TO_ANALYZE);

% --- Data Extraction and Processing Loop ---
CSP = []; % Initialize outside
CSM = []; % Initialize outside
CSPun = []; % Initialize outside

% for i = 1:length(alldata)
%      currentRow = alldata(i);
%      if isfield(currentRow, 'phase') && strcmp(currentRow.phase, PHASE_TO_ANALYZE)
%          csp_1 = currentRow.csp{1};
%          csp_2 = currentRow.csp{2};
%          csp_3 = currentRow.csp{3};
%          csp_4 = currentRow.csp{4};
%          CSP = [CSP;csp_1; csp_2;csp_3; csp_4];
% 
%          csm_1 = currentRow.csm{1};
%          csm_2 = currentRow.csm{2};
%          csm_3 = currentRow.csm{3};
%          csm_4 = currentRow.csm{4};
%          CSM = [CSM;csm_1; csm_2;csm_3; csm_4];
% 
%          cspun_1 = currentRow.cspun{1};
%          cspun_2 = currentRow.cspun{2};
%          cspun_3 = currentRow.cspun{3};
%          cspun_4 = currentRow.cspun{4};
%          CSPun = [CSPun;cspun_1; cspun_2; cspun_3; cspun_4];            
%      end
% end

for i = 1:length(alldata)
     currentRow = alldata(i);
     if isfield(currentRow, 'phase') && strcmp(currentRow.phase, PHASE_TO_ANALYZE)
        % This section now dynamically handles any number of entries in each cell array.
        % It uses vertcat to efficiently append all trials from the cell to the main matrix.
        if isfield(currentRow, 'csp') && ~isempty(currentRow.csp)
            CSP = [CSP; vertcat(currentRow.csp{:})];
        end
        
        if isfield(currentRow, 'csm') && ~isempty(currentRow.csm)
            CSM = [CSM; vertcat(currentRow.csm{:})];
        end
        
        if isfield(currentRow, 'cspun') && ~isempty(currentRow.cspun)
            CSPun = [CSPun; vertcat(currentRow.cspun{:})];
        end
     end
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
 
    time = linspace(-10, 40, size(CSP, 2));
    Rew_Early_completed = num2str(size(CSP,1));
    Rew_Early_Trial_label = ['CSP (', Rew_Early_completed,')'];


    CSM_Early_completed = num2str(size(CSM,1));
    CSM_Early_Trial_label = ['CSM (', CSM_Early_completed,')'];
  
    CSPun_Early_completed = num2str(size(CSPun,1));
    CSPun_Early_Trial_label = ['CS Pun (', CSPun_Early_completed,')'];
    
    

            %------STATS -----------------------------------------------------------------------

% %bootstrap 
tmp = bootstrap_data(CSP, 5000, 0.001);
btsrp.csp = CIadjust(tmp(1,:),tmp(2,:),tmp,size(CSP, 1),2);

tmp = bootstrap_data(CSM, 5000, 0.001);
btsrp.csm = CIadjust(tmp(1,:),tmp(2,:),tmp,size(CSM, 1),2);

tmp = bootstrap_data(CSPun, 5000, 0.001);
btsrp.cspun = CIadjust(tmp(1,:),tmp(2,:),tmp,size(CSPun, 1),2);


% %permutation tests
[perm.csp_csm, ~] = permTest_array(CSP, CSM, 1000);
[perm.csm_cspun, ~] = permTest_array(CSM, CSPun, 1000);
[perm.csp_cspun, ~] = permTest_array(CSP, CSPun, 1000);


%-----------------------------------------------------Plot CSP  -----------------------------------------------------------------------
 
%Statistical parameters
p = 0.01;
thres = 8;

figure
datasets = {CSP, CSM, CSPun};
labels = {Rew_Early_Trial_label, CSM_Early_Trial_label, CSPun_Early_Trial_label};
colors = {green, orange, red};

boottests = {'csp', 'csm', 'cspun'};
offsetsboot = [0, -0.1, -0.2, -0.3, -0.4, -0.5];
permtests = {'csp_csm', 'csm_cspun' , 'csp_cspun'};
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
line([20, 20], [yLimits(1), yLimits(2)], 'Color', 'black', 'HandleVisibility', 'off');
title(' Phase 3');  

