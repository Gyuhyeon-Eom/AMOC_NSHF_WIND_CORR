close all;
clear all;
clc;

%%
folderPath = 'C:\Users\USER\iCloudDrive\학교\4학년\1학기\PBL\Pbl'; 
%folderPath = 'C:\Users\turtl\iCloudDrive\학교\4학년\1학기\PBL\Pbl';

amocFiles = dir(fullfile(folderPath, 'AMOC*.txt'));
nshfFiles = dir(fullfile(folderPath, 'nshf*.nc'));

numModels = length(amocFiles);
amocTrends = zeros(numModels, 1);
nshfAnnualMeans = zeros(31, 121, numModels); 
legendNames = cell(1, numModels);

for i = 1:numModels
    currentFile = amocFiles(i).name;
    amocData = load(fullfile(folderPath, currentFile));
    if any(isnan(amocData(:)))
        amocData = fillmissing(amocData, 'linear');
    end
    amocAnnual = mean(reshape(amocData, 12, []), 1); 
    coeffs = polyfit(1980:2060, amocAnnual, 1); 
    amocTrends(i) = coeffs(1); 
    legendNames{i} = strrep(strrep(currentFile, 'AMOC_CMIP6_', ''), '_ssp245_1980_2060_MYM.txt', '');
end

for i = 1:numModels
    currentFile = nshfFiles(i).name;
    nshfData = ncread(fullfile(folderPath, currentFile), 'nshf');
    if any(isnan(nshfData(:)))
        nshfData = fillmissing(nshfData, 'linear', 'EndValues', 'nearest');
    end
    nshfData = reshape(permute(nshfData, [2, 3, 1, 4, 5]), 12, [], 31, 121);
    nshfAnnualMeans(:, :, i) = squeeze(mean(mean(nshfData, 1), 2)); 
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

correlations = zeros(31, 121);
pValues = zeros(31, 121);

for lat = 1:31
    for lon = 1:121
        % 모든 모델에 대해 nshf 연평균 계산
        nshfAvg = squeeze(mean(nshfAnnualMeans(lat, lon, :), 1));
        % nshfAvg와 amocTrends의 상관 계수 계산
        [R, P] = corrcoef(nshfAvg, amocTrends);
        correlations(lat, lon) = R(1, 2);
        pValues(lat, lon) = P(1, 2);
    end
end

figure;
imagesc(correlations);
colorbar;
title('Correlation between NSHF Averages and AMOC Trends Across Models');
xlabel('Longitude Index');
ylabel('Latitude Index');
axis xy;

figure;
imagesc(pValues < 0.05); % 유의 수준 0.05 이하인 그리드 셀 표시
colormap([1 1 1; 0 0 0]); % 유의미한 결과는 검은색, 그렇지 않은 경우는 흰색
title('Significant Correlations (p < 0.05)');
xlabel('Longitude Index');
ylabel('Latitude Index');
axis xy;


