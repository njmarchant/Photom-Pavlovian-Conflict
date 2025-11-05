%% Nathan Marchant July 2024
% Written for the conflict task
% A script to make changes in traces, combine them, or add extra info
% zScoring of the traces is also added at the end of this script
clear all
close all

%% define where the stuff is

tankfolder = 'C:\Photometry\PavConf\DrPhotom_Extracted\Alcohol conditioning V4 photometry_SHOCK';

%% define timestamps you wanna work with
%!!!! REMEMBER TO CHECK THE NAMES OF THE VARIABBLES ON SESDAT
var1 = {'LC+','RC+'}; %CS+ reward, all boxes
var2 = {'LC-','RC-'}; %CS-, all boxes
var3 = {'LSh','RSh'}; %CSpun, all boxes
var4 = {'LMg','RMg'}; %Magazine , all boxes



%% load individual session data
filePath = fullfile(tankfolder);
files = dir(fullfile(filePath));
files(ismember({files.name}, {'.', '..'})) = [];
files = files(~[files.isdir]);
for i = 1:length(files) %iterate through experiment folder
%     if strfind(files(i).name, r)
    if isfile(fullfile(filePath,  [files(i).name]))
        load(fullfile(filePath,  [files(i).name]))

    names = strsplit(files(i).name, {'-' , '_', ' ', '.'}); %divide the file name into separte character vectors
    r = sesdat.rat; %rat number
% Added the section below to define value in the second column of each trace 
% This is to make sure that alc=1, soc=2, trial=3, inact=4
    v1 = [];
    v1A = [];
    v1S = [];
    v1T = [];
    v1I = [];
    for j = 1:length(var1)
        if any(contains(sesdat.traces(:,1), var1(j)))
            v1A = cell2mat(sesdat.traces(contains(sesdat.traces(:,1), var1), 2));
            v1A(:,2) = 1;
        end
    end
    
    for j = 1:length(var2)
        if any(contains(sesdat.traces(:,1), var2(j)))
            v1S = cell2mat(sesdat.traces(contains(sesdat.traces(:,1), var2),2));
            v1S(:,2) = 2;
        end
    end
    
    for j = 1:length(var3)
        if any(contains(sesdat.traces(:,1), var3(j)))
            v1T = cell2mat(sesdat.traces(contains(sesdat.traces(:,1), var3),2));
            v1T(:,2) = 3;
        end
    end
    
    for j = 1:length(var4)
        if any(contains(sesdat.traces(:,1), var4(j)))
            v1I = cell2mat(sesdat.traces(contains(sesdat.traces(:,1), var4),2));
            v1I(:,2) = 4;
        end
    end
    
    arrays = {v1A,v1S,v1T,v1I};
    for j = 1:length(arrays)
        if ~isempty(arrays{j})
            if isempty(v1)
                v1= arrays{j};
            else
                v1 = cat(1, v1, arrays{j});
            end
        end
    end

       
% --- NEW SECTION: Identify CSPs that terminate in a shock ---
% This section finds each shock event (code 3) and marks the immediately
% preceding CSP event (code 1) by placing a '1' in the third column.
% This block replaces the previous loop that attempted similar logic.
if ~isempty(v1)
    % Ensure the third column exists and is initialized to 0.
    if size(v1, 2) < 3
        v1(:, 3) = 0;
    end
    
    % Get the row indices for all shock and CSP events
    shock_indices = find(v1(:, 2) == 3);
    csp_indices = find(v1(:, 2) == 1);
    
    if ~isempty(shock_indices) && ~isempty(csp_indices)
        fprintf('Found %d shocks. Checking for preceding CSPs to mark in column 3...\n', length(shock_indices));
        
        % Loop through each shock event to find its corresponding CSP
        for s_idx = 1:length(shock_indices)
            current_shock_ts = v1(shock_indices(s_idx), 1);
            
            % Find all CSPs that occurred before this shock
            preceding_csp_indices = csp_indices(v1(csp_indices, 1) < current_shock_ts);
            
            if ~isempty(preceding_csp_indices)
                % Find the index of the CSP that is closest in time before the shock
                [~, max_idx] = max(v1(preceding_csp_indices, 1));
                closest_csp_original_index = preceding_csp_indices(max_idx);
                
                % Mark this specific CSP by putting a '1' in the third column
                v1(closest_csp_original_index, 3) = 1;
                fprintf('  -> Marked CSP at index %d (time: %.2f) as shock-paired.\n', closest_csp_original_index, v1(closest_csp_original_index, 1));
            end
        end
    end
end


sesdat.traces_updated = v1;


%% 
% Here we perform z-score calulations of the traces. 
% Trace data is taken from 'sesdat.traces_updated' and after conversion
% this data is saved in 'sesdat.traces_z';

% variables to save
    sesdat.traces_z = [];
% get variables
    data = sesdat.traces_updated(:, 5:end);
    times = sesdat.traces_updated(:, 1:4);
    % dat.rat = r;
    % dat.sesdat = sesdat;
       
%define time, baseline
% Ensure that this matches what is extracted in the first script!
    time = linspace(-10, 20, size(data, 2)); %time vector the size of trace
    base = (time >= -10) & (time <= -5); %This is the time period of the trace for which the baseline is calculated
     
% zscore standardisation on df/f
    zdata = zeros(size(data));
    zbase = zeros(size(data));
    tmp = 0;
        for m = 1:size(data, 1)
            zb = mean(data(m, base));
            zsd = std(data(m,base));
            for j = 1:size(data,2)
                tmp = tmp+1;
                zbase(m, tmp) = (data(m,j) -zb);
                zdata(m,tmp) = (data(m,j) - zb)/zsd;
            end
            tmp = 0;
        end
    traces_z = [times, zdata];  
    sesdat.traces_z = [sesdat.traces_z;traces_z];

% save file (overwrites previous file, so be careful not to change
% anything)
    save([filePath '\' files(i).name(1:end-4) '.mat'], 'sesdat')

    end
end