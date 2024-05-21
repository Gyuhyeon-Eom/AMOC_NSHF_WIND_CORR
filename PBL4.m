clear all;
close all;
clc;

folderPath = 'C:\Users\turtl\iCloudDrive\학교\4학년\1학기\PBL\Pbl';
%folderPath = '/Users/eomgyuhyeon/Library/Mobile Documents/com~apple~CloudDocs/학교/4학년/1학기/PBL/Pbl';

amocFiles = dir(fullfile(folderPath, 'AMOC*.txt'));
nshfFiles = dir(fullfile(folderPath, 'nshf*.nc'));

numModels = length(amocFiles);

damoc = zeros(numModels, 1);
dnshf = zeros(121, 31, numModels);  
startYear = 1990;
endYear = 2050;
numYears = endYear - startYear + 1;

for i = 1:numModels
    currentFile = amocFiles(i).name;
    amoc_data = load(fullfile(folderPath, currentFile));
    damoc(i) = amoc_data(end-9) - amoc_data(11);

    currentFile = nshfFiles(i).name;
    nshf_data = ncread(fullfile(folderPath, currentFile), 'nshf');
    
   
    yearly_averages = zeros(size(nshf_data, 1), size(nshf_data, 2), numYears);
    for year = 1:numYears
        yearly_averages(:, :, year) = mean(nshf_data(:, :, (year-1)*12+(1:12)), 3); 
    end

    dnshf(:,:,i) = yearly_averages(:, :, end) - yearly_averages(:, :, 1);
end


correlation_matrix = zeros(31, 121);

for lat = 1:31
    for lon = 1:121
        dnshf_squeeze = squeeze(dnshf(lon, lat, :));
        correlation_matrix(lat, lon) = corr(dnshf_squeeze, damoc);
    end
end

latRange = 30:2:90;  
lonRange = -180:2:60;  

figure;  % Opens a new figure window
imagesc(lonRange, latRange, correlation_matrix);
colorbar;  % Shows the color scale
title('Correlation between DNHSF and DAMOC');  
xlabel('Longitude');  
ylabel('Latitude');  


axis tight;


colormap(jet)

daspect([1 1 1]);



%값이 아직 작다....
%연간으로 나누지 않기 
%nshf을 x축 amoc을 y축
%회귀선을 그려보고 이를 식으로 나타내자
%수식으로 검증x
%라브라도 sea의 경도 위도를 박스로 잡아서 그 지역의 평균한 nshf변화를 구하면 된다.
%모델 별로 