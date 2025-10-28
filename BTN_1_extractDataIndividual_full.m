%% Isis Alonso 18-11-20 photometry extraction script
%% --- MODIFIED: v6 - Forcing Date2 column to numeric to prevent errors ---
% This script automatically iterates through session folders, extracts photometry
% data, and uses an external Excel file to look up experimental metadata
% (phase, box) for each rat on a given date.

clear all
close all

%% ------------------- SCRIPT SETTINGS -------------------
% --- Define where the data is ---
exp = 'DrPhotom02'; 
tankfolder = 'C:\Photometry\PavConf\DrPhotom_Raw\todo';
outputfolder = 'C:\Photometry\PavConf\DrPhotom_Extracted';

% --- MODIFIED: Path to your Excel file with session metadata ---
% ❗️ IMPORTANT: Update this path to your actual file.
% The Excel file must contain these exact column headers: 'Rat', 'Box', 'Program', and 'Date2'.
phaseLookupFile = 'C:\Photometry\PavConf\Book2.xlsx';

% --- Analysis Parameters ---
time = [10, 20]; % Time window for traces: [before_event, after_event] in seconds.
res = 64; % Resampling factor
lowpass_cutoff = 3; % Low-pass cut-off in Hz
filt_steepness = .95; % Filter steepness (0.5-1)
db_atten = 90; % Decibel attenuation of filters
exc = 50; % Initial signal to exclude (in seconds)

%% ------------------- INITIALIZATION -------------------
% --- MODIFIED: Load and prepare the metadata from Excel, focusing on Date2 ---
try
    % Specify expected text format for Rat and numeric for Date2
    opts = detectImportOptions(phaseLookupFile);
    opts = setvartype(opts, 'Rat', 'string'); 
    opts = setvartype(opts, 'Date2', 'double'); % Force Date2 to be read as numeric
    opts = setvartype(opts, 'session', 'string'); % MODIFIED: Read Session column as string

    phaseTable = readtable(phaseLookupFile, opts);
    fprintf('✅ Successfully loaded and parsed the metadata Excel file.\n');
catch ME
    error('Could not read the Excel file. Please check the path and format. Error: %s', ME.message);
end

% --- Get all session folders to loop through ---
allSessionFolders = dir(tankfolder);
allSessionFolders = allSessionFolders([allSessionFolders.isdir]); % Keep only directories
allSessionFolders(ismember({allSessionFolders.name}, {'.', '..'})) = []; % Remove '.' and '..'

%% ------------------- MAIN PROCESSING LOOP -------------------
% --- Loop through each session folder ---
for s = 1:length(allSessionFolders)
    session = allSessionFolders(s).name; % The folder name is the session/date
    fprintf('\n=================================================\n');
    fprintf('Processing Session Folder: %s\n', session);
    fprintf('=================================================\n');

    % --- Individual session data extraction ---
    files = dir(fullfile(tankfolder, session));
    files(ismember({files.name}, {'.', '..'})) = [];
    
    for i = 1:length(files) % Iterate through data files in the session folder
       fprintf('Extracting data from: %s \n', files(i).name);
      
      [data, events, ts, conversion] = Sim_TDTextract(fullfile(tankfolder, session, files(i).name));
      names = strsplit(files(i).name, {'-' , '_'}); % Divide file name

    % --- SEPARATE DATA FOR LEFT AND RIGHT SETUPS ---
    for l = 1:2  % Go through the two potential rats in the filename
        if l == 1
            r = names{2}; side = 'L'; 
            raw490 = cell2mat(data(contains(data(:,1), side) & contains(data(:,1), '7'), 2)); 
            raw405 = cell2mat(data(contains(data(:,1), side) & (contains(data(:,1), '0') | contains(data(:,1), '5')), 2));
            try, ev = events(startsWith(events(:,1), side), :); catch, fprintf('No events on left side\n'); continue; end
        else
            r = names{3}; side = 'R';
            raw490 = cell2mat(data(contains(data(:,1), side) & contains(data(:,1), '7'), 2));
            raw405 = cell2mat(data(contains(data(:,1), side) & (contains(data(:,1), '0') | contains(data(:,1), '5')), 2));
            try, ev = events(startsWith(events(:,1), side), :); catch, fprintf('No events on right side\n'); continue; end
        end 
        
        % --- MODIFIED: Look up metadata using the 'Date2' column ---
        try
            % 1. Convert session folder name (e.g., '201023') to a number for direct comparison
            currentDateNum = str2double(session);
            % 2. Prepare rat ID for matching (e.g., 'R10' -> "10")
            ratID_for_lookup = string(strrep(r, 'R', ''));
            
            % 3. Find the row index by matching rat ID and the numeric Date2
            rowIndex = find(phaseTable.Rat == ratID_for_lookup & phaseTable.Date2 == currentDateNum);

            if isempty(rowIndex)
                fprintf('--> WARNING: Metadata not found for Rat %s on Date %s. Skipping this animal.\n', r, session);
                continue; % Skip to the next animal
            elseif numel(rowIndex) > 1
                 fprintf('--> WARNING: Multiple metadata entries found for Rat %s on Date %s. Using the first entry.\n', r, session);
                 rowIndex = rowIndex(1);
            end
            
            % 4. Assign metadata from the table
            type = phaseTable.Program{rowIndex}; % This will be used for the folder name
            sessionValue = phaseTable.session(rowIndex);
            sesdat.phase = type;                 % Saved in the .mat file
            sesdat.box = phaseTable.Box(rowIndex); % Saved in the .mat file
            sesdat.session = sessionValue;
            fprintf('-> Found Info: Rat: %s | Date: %s | Phase: %s | Box: %d\n', r, session, sesdat.phase, sesdat.box);
        catch ME
            fprintf('--> ERROR during metadata lookup for Rat %s, Date %s. Error: %s. Skipping.\n', r, session, ME.message);
            continue;
        end
        
        sesdat.raw490 = raw490; sesdat.raw405 = raw405;
        
        %% --- PREPROCESSING ---
        tmp = ceil(conversion*exc);
        if length(raw490) <= tmp, fprintf('--> WARNING: Data is too short. Skipping.\n'); continue; end
        filt490 = raw490(tmp:end); filt405 = raw405(tmp:end);
        ts_adjt = linspace(ts(1), ts(end), size(filt490, 1))'; 
        new_ts = ones(tmp, 1);
        
        cfFinal = controlfit(filt490, filt405);
        normDat= deltaFF(filt490,cfFinal);
        [hp_normDat, ~] = hpFilter(ts, normDat);
        lp_dFF = lpFilter(hp_normDat, conversion, lowpass_cutoff, filt_steepness, db_atten);
        lp_normDat = [new_ts; lp_dFF]; % Corrected to use lp_dFF
        
        sesdat.filt490 = filt490; sesdat.filt405 = filt405; 
        sesdat.conversion = conversion; sesdat.lp_dFF = lp_dFF; sesdat.lp_normDat = lp_normDat;
        
        % --- Plot and save for quality control ---
        a = figure('Visible', 'off');
        plot(filt490); hold on; plot(filt405); plot(1:size(lp_dFF, 1), lp_dFF);
        legend('490', '405', 'dFF'); title(sprintf('Rat %s | %s | %s', r, session, strrep(type, '_', ' ')));
        
        plotFolder = fullfile(outputfolder, 'Quality_Plots');
        if ~isfolder(plotFolder), mkdir(plotFolder); end
        saveas(a, fullfile(plotFolder, sprintf('%s_%s_%s_%s_DFF.png', exp, type, session, r)));
        close(a);
      
    %% --- MAKE TRACES ---
    if ~isempty(ev) 
        for m = 1:size(ev, 1)
            tmp = cell2mat(ev(m, 2)); tmp2 = [];
            for k = 1:size(tmp, 1)
                adjts = ceil(conversion*tmp(k, 1)); 
                try
                    signal = lp_normDat(adjts-ceil(time(1)*conversion):adjts+ceil(time(end)*conversion))'; 
                    tmp2(k, :) = signal;
                catch, fprintf('Trace timestamp out of bounds. Omitting.\n'); continue; end
            end
            if ~isempty(tmp2)
                valid_indices = any(tmp2, 2); % Find rows that were successfully populated
                tmp = tmp(valid_indices,:);
                tmp2 = tmp2(valid_indices,:);
                if ~isempty(tmp)
                    tmp(:,2) = m; tmp(:,3:4) = 0;
                    ev{m,2} = [tmp, tmp2];
                else
                    ev{m,2} = [];
                end
            else
                ev{m,2} = [];
            end
        end
        sesdat.traces = ev;
                    
        %% --- Code descriptive information into sesdat ---
        sesdat.rat = r; sesdat.exp = exp; sesdat.date = session;
        
        %% --- Save the file ---
        folderName = fullfile(outputfolder, type);
        if ~isfolder(folderName), mkdir(folderName); end
        
        saveFileName = fullfile(folderName, sprintf('%s_%s_%s_data_%s.mat', exp, type, session, r));
        save(saveFileName, 'sesdat');
        fprintf('-> Successfully saved data to: %s\n', saveFileName);
        clear sesdat;
    else 
         fprintf('No events for this rat! Skipping save.\n');
    end
    end % End of rat loop (l-loop)
    close all;
    end % End of file loop (i-loop)
end % End of session loop (s-loop)

fprintf('\n✅ All sessions processed.\n');
