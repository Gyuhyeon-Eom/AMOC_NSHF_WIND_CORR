clear all;
close all;
clc

%%
%folderPath = 'C:\Users\USER\iCloudDrive\학교\4학년\1학기\PBL\Pbl';
folderPath = 'C:\Users\turtl\iCloudDrive\학교\4학년\1학기\PBL\Pbl';

amocFiles = dir(fullfile(folderPath, 'AMOC*.txt'));
nshfFiles = dir(fullfile(folderPath, 'nshf*.nc'));

numModels = length(amocFiles);
amocTrends = zeros(numModels, 1);
legendNames = cell(1, numModels);

for i = 1:numModels
    currentFile = amocFiles(i).name;
    amocData = load(fullfile(folderPath, currentFile));
    if any(isnan(amocData(:)))
        amocData = fillmissing(amocData, 'linear');
    end
    amocAnnual = mean(reshape(amocData, 12, []), 1); 
    legendNames{i} = strrep(strrep(currentFile, 'AMOC_CMIP6_', ''), '_ssp245_1980_2060_MYM.txt', '');
end

figure;
hold on;

for i = 1:numModels
    amocData = load(fullfile(folderPath, amocFiles(i).name));
    if any(isnan(amocData(:)))
        amocData = fillmissing(amocData, 'linear');
    end
    amocMean = mean(reshape(amocData, 12, []), 1);
    plot(1980:2060, amocMean, 'DisplayName', legendNames{i});
end
legend show;
xlabel('Year');
ylabel('Mean AMOC');
title('Yearly Mean AMOC by Model');
hold off;

%%

figure;
hold on;

windowSize = 10; % in years

for i = 1:numModels
    amocData = load(fullfile(folderPath, amocFiles(i).name));
    if any(isnan(amocData(:)))
        amocData = fillmissing(amocData, 'linear');
    end
   
    amocAnnual = mean(reshape(amocData, 12, []), 1);

    amocSmoothed = movmean(amocAnnual, windowSize);

    plot(1980:2060, amocSmoothed, 'DisplayName', [legendNames{i} ' (Smoothed)']);
end

legend show;
xlabel('Year');
ylabel('10-Year Smoothed Mean AMOC');
title('10-Year Smoothed Yearly Mean AMOC by Model');
hold off;



