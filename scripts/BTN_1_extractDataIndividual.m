%% Isis Alonso 18-11-20 photometry extraction script
%% This is the main extraction and plotting script for individual sessions
%it outputs a at file per session per animal, where you can find the
%downsampled and filtered data and traces 

clear all
close all



%% define where the stuff is
exp = 'PavConf'; 
% rats = ({'R1_R7' ; 'R5_R9'; 'R6_R10'});
%the left or right set ups so make sure you know this before saving it 

tankfolder = 'C:\Photometry\PavConf\Photom_raw';
type = 'Conflict';
session= '201023';

time = [10, 40];%REMEMBER TO CHANGE THIS DEPENDING ON THEE XP!!!!!
%limits for trace making - the first one is the time before epoch, second one is after (ALWAYS KEEP BOTH POSITIVE)
% if you want plots of all sessions at the end of the extraction




%% pre-processing
res = 64; % resampling factor
lowpass_cutoff = 3; %low-pass cut-off for Hz, e.g. 2
filt_steepness = .95; %how steep the filter transition is (0.5-1, default 0.85)
db_atten = 90; %decibel attenuation of filters (default = 60)
exc = 50; %signal to take out for pre-processing (for me is 300s = 5mins)


%% individual session data extraction 
files = dir(fullfile(tankfolder,session));
files(ismember({files.name}, {'.', '..'})) = [];
for i = 1:length(files) %iterate through experiment folder
%     if strfind(files(i).name, rat)
%        if ~isfile(fullfile(tankfolder, exp, type, [ exp ' ' rat(1:2) ' ' sesdat.phase ' ' files(i).name(end-12:end-2) 'session data.mat']))...
%          ||  ~isfile(fullfile(tankfolder, exp, type, [ exp ' ' rat(5:) ' ' sesdat.phase ' ' files(i).name(end-12:end-2) 'session data.mat']))
           fprintf('Extracting %10s \n', files(i).name)
          %[data, events, ts, conversion] = Sim_TDTextract([tankfolder '\' exp '\' type '\' files(i).name]); % extract data from both set ups
          [data, events, ts, conversion] = Sim_TDTextract([tankfolder '\' session '\' files(i).name]); % extract data from both set ups
          names = strsplit(files(i).name, {'-' , '_'}); %divide the file name into separte character vectors
          %col 1 = exp name
          %col 2 = left set-up rat
          %col 3 = right set-up rat
          %col 4 = date
          %col 5 = time 

%% SEPARATE DATA TANKS LEFT AND RIGHT SET UPS
% DOUBLE CHECK THE STORE NAMES FOR THE DATA TANKS!!! One set up did weird stuff 
for l = 1:2  %go through the rats
    if l == 1
        r = names{4};
        side = 'L'; % find data from left set up
        raw490 = cell2mat(data(contains(data(:,1), side) & contains(data(:,1), '7'), 2)); 
        raw405 =cell2mat(data(contains(data(:,1), side) & contains(data(:,1), '05'), 2));
        try
        ev = events(startsWith(events(:,1), side), :); 
        catch
            fprintf('No events on left side\n')
            break
        end
    else
      r = names{5}; 
      side = 'R'; %find data from right set up
      raw490 = cell2mat(data(contains(data(:,1), side) & contains(data(:,1), '7'), 2));
      raw405 =cell2mat(data(contains(data(:,1), side) & contains(data(:,1), '05'), 2));
      try
      ev = events(startsWith(events(:,1), side), :);
      catch
        fprintf('No events on right side\n')
        break
      end
    end 

    
    %save raw traces without processing
    sesdat.raw490 = raw490;
    sesdat.raw405 = raw405;

    %% PREPROCESSING
    % CLEAN UP - exclude beginning of session when the rat is tethered/5 min wait
    tmp = ceil(conversion*exc);
    filt490 = raw490(tmp:end);
    filt405 = raw405(tmp:end);
    ts_adjt = linspace(ts(1), ts(end),  size(filt490, 1))'; %change time vector for all session based on begining exclusion
    new_ts = ts(1:tmp);
    new_ts(1:end) = ones; %this is to combine with the final data for the trace making
    clear tmp

% filt490 = raw490;
% filt405 = raw405;

    % Fit control to signal and calculate DFF
    cfFinal = controlfit(filt490, filt405);
    normDat= deltaFF(filt490,cfFinal);

%   Highpass filter and moving average...
    [hp_normDat, mov_normDat] = hpFilter(ts, normDat);

    %lowpassfilter 
    lp_normDat = lpFilter(hp_normDat, conversion, lowpass_cutoff,...
    filt_steepness, db_atten);

    lp_dFF = lp_normDat; %this is the DFF with the exclusion at the begining
    lp_normDat = [new_ts; lp_normDat]; % this is the begining as ones
    sesdat.filt490 = filt490; %save filtered 409
    sesdat.filt405 = filt405; %save filtered 405

    sesdat.conversion = conversion;
    sesdat.lp_dFF = lp_dFF; 
    sesdat.lp_normDat = lp_normDat;

%     plot and save for quality control
    a = figure;
    plot(filt490)
    hold on
    plot(filt405)
    hold on
    plot(1:size(lp_dFF), lp_dFF);
    legend('490', '405', 'dFF')
    saveas(a, [tankfolder '\' session, '\' exp ' ' type ' date ' names{6} ' DFF ' r '.png'])
  

%% MAKE TRACES...
if ~isempty(ev) %if events variable is not empty
n = 1;
for m = 1:size(ev, 1)
    tmp = cell2mat(ev(m, 2));
    for k = 1:size(tmp, 1)
        adjts = ceil(conversion*tmp(k, 1)); % get adjusted timestamps based on the sampling rate
        try
        signal = lp_normDat(adjts-ceil(time(1)*conversion):adjts+ceil(time(end)*conversion))'; 
        tmp2(k, :) = signal;
        catch
        fprintf('Trace around the timestamp goes over length of session signal! Ommiting and deleting ts...\n')
            continue %jump to the iteration
        end
    end
    if exist('tmp2', 'var') 
        if size(tmp, 1) > size(tmp2, 1)
        tmp(end-(size(tmp, 1) - size(tmp2, 1))+1:end, :) = [];
        end 
        tmp(:,2) = ones*m; tmp(:,3:4) = zeros;
    ev{m,2} = [tmp , tmp2];
    clear tmp tmp2
    end
end
sesdat.traces = ev; %save

            
%% _____________ Code descriptive information into sesdat _______________
sesdat.rat = r;
sesdat.exp = exp;
sesdat.phase = type;
sesdat.date = files(i).name(end-12:end-7);

%% - Save the file
folderName = strcat(tankfolder, '\Extracted');
    if ~isfolder(folderName)  % Check if the folder does not exist
        mkdir(folderName);  % Create the folder
        save([folderName '\' exp ' ' names{6} ' data ' r  '.mat'], 'sesdat')
    else
        save([folderName '\' exp ' ' names{6} ' data ' r  '.mat'], 'sesdat')
    end

clear sesdat 
else 
     fprintf('No events for session! Skipping to next session \n')
    continue
end
end %end of rat block
%            ses_num = ses_num +1;
% end%end of session block
           close all
end 


 

