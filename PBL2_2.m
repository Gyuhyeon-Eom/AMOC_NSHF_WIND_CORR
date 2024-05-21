close all;
clear all;
clc;

%folderPath = 'C:\Users\USER\iCloudDrive\학교\4학년\1학기\PBL\Pbl'; 
folderPath = 'C:\Users\turtl\iCloudDrive\학교\4학년\1학기\PBL\Pbl';
%folderPath = '/Users/eomgyuhyeon/Library/Mobile Documents/com~apple~CloudDocs/학교/4학년/1학기/PBL/Pbl';



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


%lat,lon 사이즈가 2x2
%month를 연평균으로 처리해도된다.
%AMOC과 NSHF의 상호작용에서 어떤게 우선순윈지 모른다.

%sensible heat flux, latent heat flux 
%shf,lhf가 양수가 되니 amoc이 약해졌다고 한방향으로 해석하는 것은 좋지 않다.
%어느지역의 shf, lhf가 amoc의 변화에 가장 민감한지를 멀티모델 분석으로 확인해보기

%어느지역의 nshf가 모형간의 amoc차이를 가장 많이 만들어 내는지를 확인해보자 
%AMOC이 강해지거나 약해졌을때 변동성이 가장 큰 지역의 HF가 가장 많이 바뀐다. 
%grid point마다 correlation을 계산한 것 처럼 grid point 마다 스무딩한 nshf 감소 y축은 amoc의 감소
%30개 모형의 자료를 비교 
%그렇게 해서 correlation map을 그려보자 
%smoothing을 안하면 급작스러운 변화가 반영이 되어서 좋은 결과를 얻기 힘들다.
%시계열을 smoothing한 후 2050-1990을 빼는게 가장 편하다. 

%grid point 마다의 x축의 nshf변화량 y축은 amoc의 변화량의 correlation을 구한뒤
%그 값을 grid point에 넣어서 correlation mapping을 해보는 것
%아마 심층수가 형성되는 지역을 잡아 낼 수 도 있다. 
%일단은 지역을 골라야 한다. 어느지역을 깊이 있게 분석해야하는지

%%
%해수 온도 자료가 있다면 이를 확인 할 수 있다.
%bulk surface evaporation equation를 가지고 수식으로 확인할 수 있다.
%SST가 1도 낮아지면 열을더 많이 흡수 할 수 있다. 
%SHF,LHF
%https://www.google.com/url?sa=i&url=https%3A%2F%2Fwww.mdpi.com%2F2072-4292%2F12%2F11%2F1796&psig=AOvVaw1wJDje647lMv47G9x2bRP9&ust=1714628228518000&source=images&cd=vfe&opi=89978449&ved=0CBIQjRxqFwoTCNiV_4Xe64UDFQAAAAAdAAAAABAE

%{
어느정도 인과관계를 corr로 증명할 수 있지만 열염순환이 약해지면 바닷물 표면 온도가 낮아짐 열순환이 줄어들면서
대기로부터 열을 더 많이 흡수 할 수도 있다. 따라서 대기가 일방적으로 바다에 주는 영향이라 해석하는 것은 위험하다.
차원을 줄였기 때문에 (2차원 자료를 1차원 자료로 줄였음) 때문에  corr가 높아진 것도 있다. 기간도 너무 길다 (60년)
대기가 해양에 주는 효과 뿐만 아니라 해양이 열염순환이 약해지며 대기의 열을 더 많이 흡수하게 되는 되먹음 작용에 의한
결과 까지도 원인으로 들어와서 해석이될 위험성이 있다.
인과관계를 확인할 수 있는 통계방법 중 하나인 lag correlation을 해봐서 1980년 부터 2040년 까지 lag
correlation을 해서 nshf가 몇년정도 앞서서 해야 correlation이 더 커지는 지를 알 수 있다.
%}