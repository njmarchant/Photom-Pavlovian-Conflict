%% Nathan Marchant May 2025
% Written for Pavlovian conflict task
% Includes: plot individual subject traces for each cue type.

% --- Initialization ---
clear all;
close all; % Good practice to close figures
dataFile = 'C:\Photometry\PavConf\IA_EtOH07AllData.mat'; % Make sure this is the correct file
load(dataFile); % This should load the 'alldata' structure

% Define the phase you want to analyze
PHASE_TO_ANALYZE = 'Conflict'; % <--- SET YOUR DESIRED PHASE HERE ('Reward', 'Conflict', 'Discrimination')

fprintf('Loading data from: %s\n', dataFile);
if ~exist('alldata', 'var')
    error('"alldata" structure not found in the loaded file. Check file content or name.');
end
fprintf('Data loaded. Processing for phase: %s\n', PHASE_TO_ANALYZE);

% --- Data Extraction and Processing per Subject ---

% Get unique subject IDs from the 'rat' field
subjects = unique({alldata.rat});
numSubjects = length(subjects);

% Initialize a struct array to hold data for each subject
% This creates a structure where each element corresponds to one subject.
subjectData = repmat(struct('rat', '', 'csp', [], 'csm', [], 'cspun', []), numSubjects, 1);

% Loop through the main 'alldata' structure to sort trials by subject
for i = 1:length(alldata)
    currentRow = alldata(i);
    
    % Check if the current row matches the desired phase
    if isfield(currentRow, 'phase') && strcmp(currentRow.phase, PHASE_TO_ANALYZE)
        % Find the index corresponding to the current subject
        subjectIdx = find(strcmp(subjects, currentRow.rat));
        
        if ~isempty(subjectIdx)
            % Store the rat ID if it's the first time we see this subject
            if isempty(subjectData(subjectIdx).rat)
                subjectData(subjectIdx).rat = subjects{subjectIdx};
            end

            % Append trial data to the correct subject's fields
            if isfield(currentRow, 'csp') && ~isempty(currentRow.csp)
                subjectData(subjectIdx).csp = [subjectData(subjectIdx).csp; vertcat(currentRow.csp{:})];
            end
            if isfield(currentRow, 'csm') && ~isempty(currentRow.csm)
                subjectData(subjectIdx).csm = [subjectData(subjectIdx).csm; vertcat(currentRow.csm{:})];
            end
            if isfield(currentRow, 'cspun') && ~isempty(currentRow.cspun)
                subjectData(subjectIdx).cspun = [subjectData(subjectIdx).cspun; vertcat(currentRow.cspun{:})];
            end
        end
    end
end

% Calculate the mean trace for each subject for each cue type
for i = 1:numSubjects
    if ~isempty(subjectData(i).csp)
        subjectData(i).csp_mean = mean(subjectData(i).csp, 1);
    else
        subjectData(i).csp_mean = []; % Handle cases where a subject has no data
    end
    if ~isempty(subjectData(i).csm)
        subjectData(i).csm_mean = mean(subjectData(i).csm, 1);
    else
        subjectData(i).csm_mean = [];
    end
    if ~isempty(subjectData(i).cspun)
        subjectData(i).cspun_mean = mean(subjectData(i).cspun, 1);
    else
        subjectData(i).cspun_mean = [];
    end
end

% --- Plotting ---

% Define time vector (assumes all data has the same number of time points)
timePoints = 0;
% Find a valid data trace to determine the number of time points
for i = 1:numSubjects
    if ~isempty(subjectData(i).csp_mean)
        timePoints = length(subjectData(i).csp_mean);
        break; % Found a valid trace, so we can exit the loop
    end
end
if timePoints == 0; error('No valid data found to create plots.'); end;
time = linspace(-10, 40, timePoints);

% Define a set of distinct colors for the plots
plotColors = lines(numSubjects);

%% Plot 1: CSP Individual Traces
figure('Name', 'CSP per Subject');
hold on;
for i = 1:numSubjects
    if ~isempty(subjectData(i).csp_mean)
        plot(time, subjectData(i).csp_mean, 'Color', plotColors(i,:), 'LineWidth', 1.5);
    end
end
hold off;
% Add reference lines for cue and outcome
ax = gca;
yLimits = ax.YLim;
line([0, 0], yLimits, 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
line([20, 20], yLimits, 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
line(ax.XLim, [0, 0], 'Color', 'k', 'LineStyle', '-');
% Add labels and title
title(['CSP Responses by Subject - ' PHASE_TO_ANALYZE ' Phase']);
xlabel('Time from Cue Onset (s)');
ylabel('dF/F (Z-score)');
legend(subjects, 'Location', 'northwest');
grid on;

%% Plot 2: CSM Individual Traces
figure('Name', 'CSM per Subject');
hold on;
for i = 1:numSubjects
    if ~isempty(subjectData(i).csm_mean)
        plot(time, subjectData(i).csm_mean, 'Color', plotColors(i,:), 'LineWidth', 1.5);
    end
end
hold off;
% Add reference lines
ax = gca;
yLimits = ax.YLim;
line([0, 0], yLimits, 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
line([20, 20], yLimits, 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
line(ax.XLim, [0, 0], 'Color', 'k', 'LineStyle', '-');
% Add labels and title
title(['CSM Responses by Subject - ' PHASE_TO_ANALYZE ' Phase']);
xlabel('Time from Cue Onset (s)');
ylabel('dF/F (Z-score)');
legend(subjects, 'Location', 'northwest');
grid on;

%% Plot 3: CSPun Individual Traces
figure('Name', 'CSPun per Subject');
hold on;
for i = 1:numSubjects
    if ~isempty(subjectData(i).cspun_mean)
        plot(time, subjectData(i).cspun_mean, 'Color', plotColors(i,:), 'LineWidth', 1.5);
    end
end
hold off;
% Add reference lines
ax = gca;
yLimits = ax.YLim;
line([0, 0], yLimits, 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
line([20, 20], yLimits, 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
line(ax.XLim, [0, 0], 'Color', 'k', 'LineStyle', '-');
% Add labels and title
title(['CSPun Responses by Subject - ' PHASE_TO_ANALYZE ' Phase']);
xlabel('Time from Cue Onset (s)');
ylabel('dF/F (Z-score)');
legend(subjects, 'Location', 'northwest');
grid on;