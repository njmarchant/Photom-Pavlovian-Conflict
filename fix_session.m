%% MATLAB Script to Update Extracted Files with Session Number
% This script iterates through existing .mat files, reads metadata from an
% Excel file, and adds a 'session' field to the 'sesdat' structure in each file.

clear all;
close all;

%% ------------------- USER SETTINGS -------------------
% --- Define the folder containing your already-extracted .mat files ---
dataFolder = 'C:\Photometry\PavConf\Conflict_Extracted\Pavlovian conditioning conflict'; % <--- SET YOUR FOLDER PATH HERE

% --- Define the path to your Excel metadata file ---
% ❗️ IMPORTANT: This file must contain the columns 'Rat', 'Date2', and 'Session'.
excelFile = 'C:\Photometry\PavConf\Book1.xlsx'; % <--- SET YOUR EXCEL FILE PATH HERE

%% ------------------- INITIALIZATION -------------------
% --- Load Metadata from Excel ---
try
    opts = detectImportOptions(excelFile);
    % Ensure columns are read with the correct data type
    opts = setvartype(opts, 'Rat', 'string');
    opts = setvartype(opts, 'Date2', 'double');
    opts = setvartype(opts, 'session', 'string'); % MODIFIED: Read Session column as string
    opts = setvartype(opts, 'session_conf', 'string'); % MODIFIED: Read Session column as string
    
    metadataTable = readtable(excelFile, opts);
    fprintf('✅ Successfully loaded metadata from: %s\n', excelFile);
catch ME
    error('Could not read the Excel file. Please check the path and column names. Error: %s', ME.message);
end

% --- Get list of all .mat files to process ---
files = dir(fullfile(dataFolder, '*.mat'));
fprintf('Found %d .mat files to process in the target folder.\n\n', length(files));

%% ------------------- MAIN PROCESSING LOOP -------------------
updated_count = 0;
skipped_count = 0;

for k = 1:length(files)
    currentFilename = files(k).name;
    fullFilePath = fullfile(dataFolder, currentFilename);
    fprintf('Processing: %s\n', currentFilename);
    
    % --- Extract Rat ID and Date from the filename ---
    % This pattern looks for a 6-digit date and a rat ID like 'R10'
    tokens = regexp(currentFilename, '_(\d{6})_data_(R\d+)\.mat$', 'tokens', 'once');
    
    if isempty(tokens)
        fprintf('  -> SKIPPING: Could not parse date and rat ID from filename.\n\n');
        skipped_count = skipped_count + 1;
        continue;
    end
    
    date_str = tokens{1};
    rat_id_str = tokens{2};
    
    % --- Prepare data for matching with the Excel table ---
    date_num = str2double(date_str);
    rat_id_for_lookup = string(strrep(rat_id_str, 'R', '')); % 'R10' -> "10"
    
    % --- Find the matching row in the metadata table ---
    rowIndex = find(metadataTable.Date2 == date_num & metadataTable.Rat == rat_id_for_lookup);
    
    if isempty(rowIndex)
        fprintf('  -> SKIPPING: No matching entry found in Excel for Rat %s on Date %s.\n\n', rat_id_str, date_str);
        skipped_count = skipped_count + 1;
        continue;
    elseif numel(rowIndex) > 1
        fprintf('  -> WARNING: Multiple entries found. Using the first one.\n');
        rowIndex = rowIndex(1);
    end
    
    % --- Get the session value from the table ---
    sessionValue = metadataTable.session(rowIndex);
    session_C_Value = metadataTable.session_conf(rowIndex);
    
    % --- Load, Update, and Save the .mat file ---
    try
        load(fullFilePath); % Loads the 'sesdat' structure
        
        if exist('sesdat', 'var')
            % Add the new field
            sesdat.session = sessionValue;
            sesdat.session_conf = session_C_Value;
            
            % Save the updated structure back to the file
            save(fullFilePath, 'sesdat');
            % MODIFIED: Removed num2str to correctly print the string value
            fprintf('  -> SUCCESS: Added session ''%s'' to the file.\n\n', sessionValue);
            updated_count = updated_count + 1;
        else
            fprintf('  -> SKIPPING: File does not contain a "sesdat" variable.\n\n');
            skipped_count = skipped_count + 1;
        end
        
    catch ME
        fprintf('  -> ERROR: Could not load or save the file. Error: %s\n\n', ME.message);
        skipped_count = skipped_count + 1;
    end
    
    % Clear the loaded variable to prevent conflicts
    clear sesdat;
end

fprintf('\n--- Update Complete ---\n');
fprintf('Successfully updated %d files.\n', updated_count);
fprintf('Skipped %d files.\n', skipped_count);
