%% Nathan Marchant August 2025
% This script calculates the mean z-scored photometry signal in 1-second time
% bins for each subject and session from a folder of data files.
% It saves detailed data per rat, including pooled early/late session averages,
% and a final collated summary file.

close all;
clear all;

%% --- 1. Setup ---
% Define the folder containing your individual session .mat files
dataFolder = 'C:\Photometry\PavConf\DrPhotom_Extracted\Alcohol conditioning V4 photometry';

% --- MODIFIED: Define which sessions to group together ---
% List the session names you want to pool for the "first" group
first_sessions = {'retrain_1', 'retrain_2'};
% List the session names for the "final" group for most rats
final_sessions_default = {'retrain_4', 'retrain_5'};
% List the session names for the "final" group for special case rats
final_sessions_special = {'retrain_4'};
special_rats = {'R15', 'R16'}; % List of rats that use the special final group

% --- Internal Setup ---
ratDataMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
collatedResultsCell = {};
fprintf('Starting data collation from folder: %s\n', dataFolder);

%% --- 2. Data Collation Loop ---
files = dir(fullfile(dataFolder, '*.mat'));

for i = 1:length(files)
    if contains(files(i).name, '_binned_means.mat') || contains(files(i).name, 'ALL_RATS_collated_means')
        continue;
    end
    
    fprintf('Loading and processing file: %s\n', files(i).name);
    load(fullfile(dataFolder, files(i).name), 'sesdat');
    
    if exist('sesdat', 'var') && all(isfield(sesdat, {'rat', 'session', 'traces_z'}))
        ratID = sesdat.rat;
        if isnumeric(ratID); ratID = ['Rat' num2str(ratID)]; end
        sessionID = sesdat.session;
        
        if ~isKey(ratDataMap, ratID)
            ratDataMap(ratID) = struct(...
                'sessionMap', containers.Map('KeyType', 'char', 'ValueType', 'any'), ...
                'early_pool', struct('csp', [], 'csm', []), ...
                'late_pool',  struct('csp', [], 'csm', []) ...
            );
        end
        ratStruct = ratDataMap(ratID);
        
        % --- Handle Individual Session Data ---
        if ~isKey(ratStruct.sessionMap, sessionID)
            ratStruct.sessionMap(sessionID) = struct('csp', [], 'csm', []);
        end
        currentSessionData = ratStruct.sessionMap(sessionID);
        
        csp_traces = sesdat.traces_z(sesdat.traces_z(:, 2) == 1, 5:end);
        csm_traces = sesdat.traces_z(sesdat.traces_z(:, 2) == 2, 5:end);
        
        
        currentSessionData.csp = [currentSessionData.csp; csp_traces];
        currentSessionData.csm = [currentSessionData.csm; csm_traces];

        ratStruct.sessionMap(sessionID) = currentSessionData;
        
        % --- Handle Pooled Session Data ---
        if ismember(sessionID, first_sessions)
            ratStruct.early_pool.csp = [ratStruct.early_pool.csp; csp_traces];
            ratStruct.early_pool.csm = [ratStruct.early_pool.csm; csm_traces];

        end
        
        applicable_late_sessions = final_sessions_default;
        if ismember(ratID, special_rats)
            applicable_late_sessions = final_sessions_special;
        end
        if ismember(sessionID, applicable_late_sessions)
            ratStruct.late_pool.csp = [ratStruct.late_pool.csp; csp_traces];
            ratStruct.late_pool.csm = [ratStruct.late_pool.csm; csm_traces];
   
        end
        
        ratDataMap(ratID) = ratStruct;
    else
        fprintf('  -> Skipped file. "sesdat" structure not found or invalid.\n');
    end
    clear sesdat;
end

fprintf('\nData collation complete. Found data for %d rats.\n\n', length(keys(ratDataMap)));

%% --- 3. Analysis and Saving Loop ---
allRatIDs = keys(ratDataMap);
time = linspace(-10, 20, 501);

for r = 1:length(allRatIDs)
    ratID = allRatIDs{r};
    ratStruct = ratDataMap(ratID);
    fprintf('Analyzing data for rat: %s\n', ratID);
    
    ratResults = struct();
    
    % --- Analyze Individual Sessions ---
    allSessionIDs = keys(ratStruct.sessionMap);
    for s = 1:length(allSessionIDs)
        sessionID = allSessionIDs{s};
        sessionData = ratStruct.sessionMap(sessionID);
        [sessionResults, collatedResultsCell] = analyzeAndCollate(sessionData, time, ratID, sessionID, collatedResultsCell);
        if ~isempty(fieldnames(sessionResults))
            ratResults.(sessionID) = sessionResults;
        end
    end
    
    % --- Analyze Pooled Early/Late Sessions ---
    [earlyResults, collatedResultsCell] = analyzeAndCollate(ratStruct.early_pool, time, ratID, 'EARLY_POOL', collatedResultsCell);
    if ~isempty(fieldnames(earlyResults)); ratResults.EARLY_POOL = earlyResults; end
    
    [lateResults, collatedResultsCell] = analyzeAndCollate(ratStruct.late_pool, time, ratID, 'LATE_POOL', collatedResultsCell);
    if ~isempty(fieldnames(lateResults)); ratResults.LATE_POOL = lateResults; end

    % Save the detailed results for the current rat
    if ~isempty(fieldnames(ratResults))
        outputFilename = fullfile(dataFolder, sprintf('%s_binned_means.mat', ratID));
        save(outputFilename, 'ratResults');
        fprintf('  -> Saved detailed results (including pools) to: %s\n', outputFilename);
    end
end

%% --- 4. Final Collation and Saving ---
if ~isempty(collatedResultsCell)
    fprintf('\nCreating final summary file...\n');
    binNames = "Bin" + (1:10);
    variableNames = [{'Rat', 'Session', 'TrialType'}, binNames];
    collatedMeans = cell2table(collatedResultsCell, 'VariableNames', variableNames);
    
    matFilename = fullfile(dataFolder, 'ALL_RATS_collated_means.mat');
    csvFilename = fullfile(dataFolder, 'ALL_RATS_collated_means.csv');
    save(matFilename, 'collatedMeans');
    writetable(collatedMeans, csvFilename);
    fprintf('Successfully saved summary file to:\n  - %s\n  - %s\n', matFilename, csvFilename);
else
    fprintf('\nNo data was processed, final summary file was not created.\n');
end

fprintf('\nAnalysis finished.\n');

%% --- Helper Function ---
function [results, collatedCell] = analyzeAndCollate(dataStruct, time, ratID, label, collatedCell)
    % This function calculates binned means and appends to the master list
    results = struct();
    trialTypes = {'csp', 'csm'};
    
    for tt = 1:length(trialTypes)
        trialType = trialTypes{tt};
        traces = dataStruct.(trialType);
        
        if isempty(traces); continue; end
        
        num_trials = size(traces, 1);
        num_bins = 10;
        binned_trial_means = NaN(num_trials, num_bins);
        
        for bin = 1:num_bins
            startTime = bin - 1;
            endTime = bin;
            time_indices = (time >= startTime) & (time < endTime);
            binned_trial_means(:, bin) = mean(traces(:, time_indices), 2);
        end
        
        mean_per_bin = mean(binned_trial_means, 1, 'omitnan');
        sem_per_bin = std(binned_trial_means, 0, 1, 'omitnan') / sqrt(num_trials);
        
        results.([trialType '_mean_per_bin']) = mean_per_bin;
        results.([trialType '_sem_per_bin']) = sem_per_bin;
        
        newRow = [{ratID, label, trialType}, num2cell(mean_per_bin)];
        collatedCell = [collatedCell; newRow];
    end
end