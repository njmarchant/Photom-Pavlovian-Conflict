%% Collate zmean scores of multiple sessions into a single file per rat 

fclose all;
clear all;
close all;

%% define where the stuff is
tankfolder = 'C:\Photometry\csmial_045 Photom\SA_Extracted_Data_230930_(-5_+10)\240924_zMean\';
phase = 'SA';

%% find files
filePath = fullfile(tankfolder,phase);
filesAndFolders = dir(fullfile(filePath));
files = filesAndFolders(~[filesAndFolders.isdir]); 
files(ismember({files.name}, {'.', '..'})) = [];
numfiles = length(files);

%% variables to save

collated.all.csp_labels = [];
collated.all.csp = [];
collated.all.csm_labels = [];
collated.all.csm = [];
collated.all.cspun_labels = [];
collated.all.cspun = [];


collated.mean.labels = cell(18,2);
collated.mean.csp = zeros(18,4);
collated.mean.csm = zeros(18,4);
collated.mean.cspun = zeros(18,4);

% Loop through files to collate data
load(fullfile(filePath, files(i).name));
for i = 1:18
    
    collated.all.csp_labels = [collated.all.csp_labels; [alldata(i).rat,alldata(i).phase]];
    collated.all.csp = [collated.all.csp; alldat.zmean_data.csp];
    
    collated.all.csm_labels = [collated.all.csm_labels; alldat.zmean_labels.csm_labels];
    collated.all.csm = [collated.all.csm; alldat.zmean_data.csm];
    

    

    if ~isempty(alldat.zmean_labels.csm_labels)
        collated.mean.labels{i,1} = alldat.zmean_labels.csm_labels{1,1};
        collated.mean.labels{i,2} = alldat.zmean_labels.csm_labels{1,2};
        collated.mean.labels{i,3} = phase;
        
    end
    % Check if the cspohol means are NaN and if not add it to the array
    if ~isnan(alldat.zmean_data.cspmean) 
         collated.mean.csp(i, :) = alldat.zmean_data.cspmean;
    end
      
    % Check if the csmial means are NaN and if not add it to the array
    if ~isnan(alldat.zmean_data.csmmean) 
         collated.mean.csm(i, :) = alldat.zmean_data.csmmean;
    end
    
    
end
   p = char(phase);
    folderName = strcat(tankfolder,phase);
        save([folderName '\csmial 045 ' p ' zMean all rats.mat'], 'collated')