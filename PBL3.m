clear all;
clc;

%% 
folderPath = 'C:\Users\turtl\iCloudDrive\학교\4학년\1학기\PBL\Pbl';

amocFiles = dir(fullfile(folderPath, 'AMOC*.txt'));
nshfFiles = dir(fullfile(folderPath, 'nshf*.nc'));

numModels = length(amocFiles);

%% 
damoc = zeros(numModels, 1); % damoc 초기화
for i = 1:numModels
    currentFile = amocFiles(i).name;
    amoc_data = load(fullfile(folderPath, currentFile));
    amoc_smoothed = smooth(amoc_data, 10);
    damoc(i) = mean(amoc_smoothed); 
end

%%
dnshf = zeros(31, 121, numModels); % dnshf 초기화
for i = 1:numModels
    currentFile = nshfFiles(i).name;
    nshf_data = ncread(fullfile(folderPath, currentFile), 'nshf');
    for lat = 1:size(nshf_data, 1)
        for lon = 1:size(nshf_data, 2)
            temp_data = squeeze(nshf_data(lat, lon, :, :));
            temp_smoothed = smooth(temp_data, 10);
            dnshf(lat, lon, i) = mean(temp_smoothed);
        end
    end
end

%%
correlation_matrix = zeros(31, 121);
for lat = 1:31
    for lon = 1:121

        correlation_matrix(lat, lon) = corr(damoc, squeeze(dnshf(lat, lon, :)));
    end
end

figure;
imagesc(correlation_matrix);
colorbar;
title('Correlation between DAMOC and DNSHF across all models');
xlabel('Longitude index');
ylabel('Latitude index');
