%% Nathan Marchant July 2024
% Written for the conflict task
% This script will calculate the zMean for a given period, defined below

close all
clear all
close all

% define where the stuff is
tankfolder = 'C:\Photometry\PavConf';

    temp.zmean = [];
    % collated_zMean = [];
    collated_zMean.collated_csp_labels = [];
    collated_zMean.collated_csp = [];
    collated_zMean.collated_csm_labels = [];
    collated_zMean.collated_csm = [];
    collated_zMean.collated_cspun_labels = [];
    collated_zMean.collated_cspun = [];

    collated_zMean.csp_Mean_labels = [];
    collated_zMean.csp_Mean = [];
    collated_zMean.csm_Mean_labels = [];
    collated_zMean.csm_Mean = [];
    collated_zMean.cspun_Mean_labels = [];
    collated_zMean.cspun_Mean = [];

% load invidivual session data
filePath = fullfile(tankfolder);
filesAndFolders = dir(fullfile(filePath));
files = filesAndFolders(~[filesAndFolders.isdir]); 
files(ismember({files.name}, {'.', '..'})) = [];
for i = 1:18
    load(fullfile(filePath,  'IA_EtOH07AllData.mat'))

    % get variables
    r = alldata(i).rat;
    phase = alldata(i).phase;
    cspdata = alldata(i).csp{1,1}(:,1:end);
    %define time, baseline etc

    time = linspace(-10, 40, size(cspdata, 2)); %time vector the size of trace
    % base = (time >= -10) & (time <= -5); %This is the time period of the trace for which the baseline is calculated
    lims1 = (time >= -5) & (time <= 0); %baseline - pre cue
    lims2 = (time >= 0) & (time <= 2); %limits - cue to alcohol
    lims3 = (time >= 15) & (time <= 20); %limits - alcohol to cue off
    lims4 = (time >= 20) & (time <= 22); %limits - shock period (2 sec ??)

    for j = 1:numel(alldata(i).csp)
        %calculate mean in CS Plus
        cspdata = alldata(i).csp{1,j}(:,1:end);
        temp.zmean.csp_pre = mean(cspdata(:, lims1), 2);
        temp.zmean.csp_cue = mean(cspdata(:, lims2), 2);
        temp.zmean.csp_alc = mean(cspdata(:, lims3), 2);
        temp.zmean.csp_shock = mean(cspdata(:, lims4), 2);
        % add rat, sec, hemi, to the respective values
        csp_r_col = [repmat({r}, size(cspdata, 1),1)];
        csp_p_col = [repmat({phase}, size(cspdata, 1),1)];
        ses_mean_pre = mean(temp.zmean.csp_pre,1);
        ses_mean_cue = mean(temp.zmean.csp_cue,1);
        ses_mean_alc = mean(temp.zmean.csp_alc,1);
        ses_mean_shock = mean(temp.zmean.csp_shock,1);

        collated_csp_labels = [csp_r_col, csp_p_col];
        collated_csp = [temp.zmean.csp_pre, temp.zmean.csp_cue, temp.zmean.csp_alc, temp.zmean.csp_shock];
        collated_zMean.collated_csp_labels = [collated_zMean.collated_csp_labels;collated_csp_labels];
        collated_zMean.collated_csp = [collated_zMean.collated_csp;collated_csp];
        
        collated_csp_Mean_labels = [csp_r_col(1,1), csp_p_col(1,1)];
        collated_csp_Mean = [ses_mean_pre, ses_mean_cue, ses_mean_alc, ses_mean_shock];
        collated_zMean.csp_Mean_labels = [collated_zMean.csp_Mean_labels;collated_csp_Mean_labels];
        collated_zMean.csp_Mean = [collated_zMean.csp_Mean;collated_csp_Mean];

    end
    for j = 1:numel(alldata(i).csm)
        %calculate mean in CS Plus
        csmdata = alldata(i).csm{1,j}(:,1:end);
        temp.zmean.csm_pre = mean(csmdata(:, lims1), 2);
        temp.zmean.csm_cue = mean(csmdata(:, lims2), 2);
        temp.zmean.csm_alc = mean(csmdata(:, lims3), 2);
        temp.zmean.csm_shock = mean(csmdata(:, lims4), 2);
        % add rat, sec, hemi, to the respective values
        csm_r_col = [repmat({r}, size(csmdata, 1),1)];
        csm_p_col = [repmat({phase}, size(csmdata, 1),1)];
        ses_mean_pre = mean(temp.zmean.csm_pre,1);
        ses_mean_cue = mean(temp.zmean.csm_cue,1);
        ses_mean_alc = mean(temp.zmean.csm_alc,1);
        ses_mean_shock = mean(temp.zmean.csm_shock,1);

        collated_csm_labels = [csm_r_col, csm_p_col];
        collated_csm = [temp.zmean.csm_pre, temp.zmean.csm_cue, temp.zmean.csm_alc, temp.zmean.csm_shock];
        collated_zMean.collated_csm_labels = [collated_zMean.collated_csm_labels;collated_csm_labels];
        collated_zMean.collated_csm = [collated_zMean.collated_csm;collated_csm];

        collated_csm_Mean_labels = [csp_r_col(1,1), csp_p_col(1,1)];
        collated_csm_Mean = [ses_mean_pre, ses_mean_cue, ses_mean_alc, ses_mean_shock];
        collated_zMean.csm_Mean_labels = [collated_zMean.csm_Mean_labels;collated_csm_Mean_labels];
        collated_zMean.csm_Mean = [collated_zMean.csm_Mean;collated_csm_Mean];

    end
    for j = 1:numel(alldata(i).cspun)
        %calculate mean in CS Plus
        cspundata = alldata(i).cspun{1,j}(:,1:end);
        temp.zmean.cspun_pre = mean(cspundata(:, lims1), 2);
        temp.zmean.cspun_cue = mean(cspundata(:, lims2), 2);
        temp.zmean.cspun_alc = mean(cspundata(:, lims3), 2);
        temp.zmean.cspun_shock = mean(cspundata(:, lims4), 2);
        % add rat, sec, hemi, to the respective values
        cspun_r_col = [repmat({r}, size(cspundata, 1),1)];
        cspun_p_col = [repmat({phase}, size(cspundata, 1),1)];
        ses_mean_pre = mean(temp.zmean.cspun_pre,1);
        ses_mean_cue = mean(temp.zmean.cspun_cue,1);
        ses_mean_alc = mean(temp.zmean.cspun_alc,1);
        ses_mean_shock = mean(temp.zmean.cspun_shock,1);

        collated_cspun_labels = [cspun_r_col, cspun_p_col];
        collated_cspun = [temp.zmean.cspun_pre, temp.zmean.cspun_cue, temp.zmean.cspun_alc, temp.zmean.cspun_shock];
        collated_zMean.collated_cspun_labels = [collated_zMean.collated_cspun_labels;collated_cspun_labels];
        collated_zMean.collated_cspun = [collated_zMean.collated_cspun;collated_cspun];

        collated_cspun_Mean_labels = [csp_r_col(1,1), csp_p_col(1,1)];
        collated_cspun_Mean = [ses_mean_pre, ses_mean_cue, ses_mean_alc, ses_mean_shock];
        collated_zMean.cspun_Mean_labels = [collated_zMean.cspun_Mean_labels;collated_cspun_Mean_labels];
        collated_zMean.cspun_Mean = [collated_zMean.cspun_Mean;collated_cspun_Mean];
        
    end
    
    
end
                


save([filePath '\IA_ETOH07_collated.mat'], 'collated_zMean')