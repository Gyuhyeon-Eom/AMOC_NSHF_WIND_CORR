clear all;
close all;
clc;

%folderPath ='/Users/eomgyuhyeon/Library/Mobile Documents/com~apple~CloudDocs/학교/4학년/1학기/PBL/Pbl';
folderPath = 'C:\Users\turtl\iCloudDrive\학교\4학년\1학기\PBL\Pbl';
amocFiles = dir(fullfile(folderPath, 'AMOC*.txt'));
nshfFiles = dir(fullfile(folderPath, 'nshf*.nc'));
windFiles = dir(fullfile(folderPath, 'sfcWind*.nc'));

numModels = 30;
startYear = 1990;
endYear = 2050;
numYears = endYear - startYear + 1;

damoc = zeros(numModels, 1);
dnshf = zeros(121, 31, numModels);   

for i = 1:numModels
    currentFile = amocFiles(i).name;
    amoc_data = load(fullfile(folderPath, currentFile));
    damoc(i) = amoc_data(end-9) - amoc_data(11);  
end

for i = 1:numModels
    currentFile = nshfFiles(i).name;
    nshfPath = fullfile(folderPath, currentFile);

    nshfData = ncread(nshfPath, 'nshf');

    startYearIndex = (startYear - 1980) * 12 + 1;  
    endYearIndex = (endYear - 1980) * 12;         


    startNshf = mean(nshfData(:, :, startYearIndex:startYearIndex+11), 3);  
    endNshf = mean(nshfData(:, :, endYearIndex-11:endYearIndex), 3);        

    dnshf(:, :, i) = endNshf - startNshf; 
end

correlation_matrix = zeros(31, 121);
for lat = 1:31
    for lon = 1:121
        dnshf_squeeze = squeeze(dnshf(lon, lat, :));
        if any(~isnan(dnshf_squeeze))
            correlation_matrix(lat, lon) = corr(dnshf_squeeze, damoc, 'Rows', 'complete');
        else
            correlation_matrix(lat, lon) = NaN;
        end
    end
end

latRange = 30:2:90;
lonRange = -180:2:60;

figure;
hold on;
imagesc(lonRange, latRange, correlation_matrix, 'AlphaData', ~isnan(correlation_matrix));
colorbar;
title('Correlation between DNHSF and DAMOC');
xlabel('Longitude');
ylabel('Latitude');
hold off;

%%
correlation_matrix = zeros(31, 121);
high_corr_indices = [];

for lat = 1:31
    for lon = 1:121
        dnshf_squeeze = squeeze(dnshf(lon, lat, :));
        if any(~isnan(dnshf_squeeze))
            corr_value = corr(dnshf_squeeze, damoc, 'Rows', 'complete');
            correlation_matrix(lat, lon) = corr_value;
            if abs(corr_value) >= 0.6
                high_corr_indices = [high_corr_indices; [lon, lat]];
            end
        else
            correlation_matrix(lat, lon) = NaN;
        end
    end
end

avg_nshf_high_corr = zeros(numModels, 1);

for i = 1:numModels
    nshf_values = [];
    for idx = 1:size(high_corr_indices, 1)
        lon = high_corr_indices(idx, 1);
        lat = high_corr_indices(idx, 2);
        nshf_values = [nshf_values, dnshf(lon, lat, i)];
    end
    avg_nshf_high_corr(i) = mean(nshf_values);
end

figure;
scatter(avg_nshf_high_corr, damoc, 'filled');
hold on;  

p = polyfit(avg_nshf_high_corr, damoc, 1); 
x_fit = linspace(min(avg_nshf_high_corr), max(avg_nshf_high_corr), 100);  
y_fit = polyval(p, x_fit);  

plot(x_fit, y_fit, 'r-', 'LineWidth', 2);  
hold off;

title('Correlation between DNHSF and DAMOC in High Correlation Areas');
xlabel('Average DNHSF Change in High Correlation Areas');
ylabel('AMOC Change');
grid on;

fprintf('The regression line is y = %.2fx + %.2f\n', p(1), p(2));

slope = p(1);
intercept = p(2);

disp(['Regression Line: y = ' num2str(slope) 'x + ' num2str(intercept)]);

%%


numModels = numel(windFiles); 
startYear = 1990;
endYear = 2050;

dnshf = zeros(121, 31, numModels);
dWind = zeros(121, 31, numModels);

for i = 1:numModels
    nshfPath = fullfile(folderPath, nshfFiles(i).name);
    nshfData = ncread(nshfPath, 'nshf');

    startYearIndex = (startYear - 1980) * 12 + 1;
    endYearIndex = (endYear - 1980) * 12;

    startNshf = mean(nshfData(:, :, startYearIndex:startYearIndex+11), 3);
    endNshf = mean(nshfData(:, :, endYearIndex-11:endYearIndex), 3);

    dnshf(:, :, i) = endNshf - startNshf;

    windPath = fullfile(folderPath, windFiles(i).name);
    windData = ncread(windPath, 'sfcWind');

    startWind = mean(windData(:, :, startYearIndex:startYearIndex+11), 3);
    endWind = mean(windData(:, :, endYearIndex-11:endYearIndex), 3);

    dWind(:, :, i) = endWind - startWind;
end

correlation_matrix_wind_nshf = zeros(31, 121);
for lat = 1:31
    for lon = 1:121
        wind_changes = squeeze(dWind(lon, lat, :));
        nshf_changes = squeeze(dnshf(lon, lat, :));

        validIdx = ~isnan(wind_changes) & ~isnan(nshf_changes);
        if any(validIdx)
            correlation_matrix_wind_nshf(lat, lon) = corr(wind_changes(validIdx), nshf_changes(validIdx));
        else
            correlation_matrix_wind_nshf(lat, lon) = NaN;
        end
    end
end

latRange = 30:2:90;
lonRange = -180:2:60;

figure;
hold on;
imagesc(lonRange, latRange, correlation_matrix_wind_nshf, 'AlphaData', ~isnan(correlation_matrix_wind_nshf));
colorbar;
title('Correlation between Wind Speed Change and NSHF Change (1990-2050)');
xlabel('Longitude');
ylabel('Latitude');
hold off;

%% wind annualmean

years = 1990:2050;
annualMeans = zeros(length(years), length(windFiles));

for i = 1:length(windFiles)
    filePath = fullfile(windFiles(i).folder, windFiles(i).name);  % Full path to the current file
    windData = ncread(filePath, 'sfcWind');

    monthlyMean = squeeze(mean(windData, [1, 2, 5], 'omitnan'));

    annualMean = mean(reshape(monthlyMean, 12, []), 1, 'omitnan');

    startYear = 1990 - ncread(filePath, 'year', 1, 1) + 1;
    endYear = 2050 - ncread(filePath, 'year', 1, 1) + 1;
    annualMeans(:, i) = annualMean(startYear:endYear);
end

figure;
hold on;
colors = lines(length(windFiles)); 

for i = 1:length(windFiles)
    plot(years, annualMeans(:, i), 'Color', colors(i, :), 'LineWidth', 2);
end
hold off;

legend({windFiles.name}, 'Interpreter', 'none', 'Location', 'best');
title('Annual Mean Surface Wind Speeds from NC Files (1990-2050)');
xlabel('Year');
ylabel('Mean Wind Speed (m/s)');

%% wind annualmean smoothing 10years

years = 1990:2050;
annualMeans = zeros(length(years), length(windFiles));

for i = 1:length(windFiles)
    filePath = fullfile(windFiles(i).folder, windFiles(i).name);  % Full path to the current file
    windData = ncread(filePath, 'sfcWind');

    monthlyMean = squeeze(mean(windData, [1, 2, 5], 'omitnan'));

    annualMean = mean(reshape(monthlyMean, 12, []), 1, 'omitnan');

    startYear = 1990 - ncread(filePath, 'year', 1, 1) + 1;
    endYear = 2050 - ncread(filePath, 'year', 1, 1) + 1;
    annualMeans(:, i) = annualMean(startYear:endYear);
end

smoothingPeriod = 10;
smoothedAnnualMeans = movmean(annualMeans, smoothingPeriod, 1);

figure;
hold on;
colors = lines(length(windFiles));  % Generate distinct colors for each model

for i = 1:length(windFiles)
    plot(years, smoothedAnnualMeans(:, i), 'Color', colors(i, :), 'LineWidth', 2);
end
hold off;

legend({windFiles.name}, 'Interpreter', 'none', 'Location', 'best');
title('Smoothed Annual Mean Surface Wind Speeds from NC Files (1990-2050)');
xlabel('Year');
ylabel('Mean Wind Speed (m/s)');

%% wind all annualmean


years = 1990:2050;

allModelsAnnualMeans = zeros(length(years), length(windFiles));

for i = 1:length(windFiles)
    filePath = fullfile(windFiles(i).folder, windFiles(i).name);  % Full path to the current file
    windData = ncread(filePath, 'sfcWind');

    monthlyMean = squeeze(mean(windData, [1, 2, 5], 'omitnan'));

    annualMean = mean(reshape(monthlyMean, 12, []), 1, 'omitnan');

    startYear = 1990 - ncread(filePath, 'year', 1, 1) + 1;
    endYear = 2050 - ncread(filePath, 'year', 1, 1) + 1;
    allModelsAnnualMeans(:, i) = annualMean(startYear:endYear);
end

meanAcrossModels = mean(allModelsAnnualMeans, 2);

figure;
plot(years, meanAcrossModels, 'b-', 'LineWidth', 2);
title('Average Annual Mean Surface Wind Speeds Across Models (1990-2050)');
xlabel('Year');
ylabel('Average Mean Wind Speed (m/s)');
grid on;

%% nshf annualmean(smoothing 10years)

years = 1990:2050;

annualMeans = zeros(length(years), length(nshfFiles));

for i = 1:length(nshfFiles)
    filePath = fullfile(nshfFiles(i).folder, nshfFiles(i).name);  % Full path to the current file
    nshfData = ncread(filePath, 'nshf');  % Read the nshf variable

    monthlyMean = squeeze(mean(nshfData, [1, 2, 5], 'omitnan'));

    annualMean = mean(reshape(monthlyMean, 12, []), 1, 'omitnan');

    startYearIndex = find(years == 1990);
    endYearIndex = find(years == 2050);
    actualStartYear = ncread(filePath, 'year', 1, 1);  % Reading the start year from the file
    startYear = 1990 - actualStartYear + 1;
    endYear = 2050 - actualStartYear + 1;
    annualMeans(:, i) = annualMean(startYear:endYear);
end

smoothingPeriod = 10;
smoothedAnnualMeans = movmean(annualMeans, smoothingPeriod, 1);

figure;
hold on;
colors = lines(length(nshfFiles));  % Generate distinct colors for each model

for i = 1:length(nshfFiles)
    plot(years, smoothedAnnualMeans(:, i), 'Color', colors(i, :), 'LineWidth', 2);
end
hold off;

legend({nshfFiles.name}, 'Interpreter', 'none', 'Location', 'best');
title('Smoothed Annual Mean Net Surface Heat Flux from NC Files (1990-2050)');
xlabel('Year');
ylabel('Mean Net Surface Heat Flux (W/m^2)');




%% nshf all annualmean

years = 1990:2050;

allModelsAnnualMeans = zeros(length(years), length(nshfFiles));

for i = 1:length(nshfFiles)
    filePath = fullfile(nshfFiles(i).folder, nshfFiles(i).name);  % Full path to the current file
    nshfData = ncread(filePath, 'nshf');

    % Calculate monthly means, omitting NaNs, and squeeze to remove singleton dimensions
    monthlyMean = squeeze(mean(nshfData, [1, 2, 5], 'omitnan'));

    % Compute annual mean for each year, assuming data covers January to December each year
    annualMean = mean(reshape(monthlyMean, 12, []), 1, 'omitnan');

    % Subset the annualMean array to include only the years 1990 to 2050
    % Adjusting to start at the first year provided by the netCDF year variable
    startYear = 1990 - ncread(filePath, 'year', 1, 1) + 1;
    endYear = 2050 - ncread(filePath, 'year', 1, 1) + 1;
    allModelsAnnualMeans(:, i) = annualMean(startYear:endYear);
end

meanAcrossModels = mean(allModelsAnnualMeans, 2);

figure;
plot(years, meanAcrossModels, 'b-', 'LineWidth', 2);
title('Average Annual Mean Net Surface Heat Flux Across Models (1990-2050)');
xlabel('Year');
ylabel('Average Mean Net Surface Heat Flux (W/m^2)');
grid on;

%% nshf, wind data
annualWindMeans = zeros(length(1990:2050), numModels);

for i = 1:numModels
    filePath = fullfile(windFiles(i).folder, windFiles(i).name);
    windData = ncread(filePath, 'sfcWind');
    monthlyMean = squeeze(mean(windData, [1, 2, 5], 'omitnan'));
    annualMean = mean(reshape(monthlyMean, 12, []), 1, 'omitnan');
    startYear = 1990 - ncread(filePath, 'year', 1, 1) + 1;
    endYear = 2050 - ncread(filePath, 'year', 1, 1) + 1;
    annualWindMeans(:, i) = annualMean(startYear:endYear);
end

nshfFiles = dir(fullfile(folderPath, 'nshf*.nc'));
annualNshfMeans = zeros(length(1990:2050), numModels);

for i = 1:numModels
    filePath = fullfile(nshfFiles(i).folder, nshfFiles(i).name);
    nshfData = ncread(filePath, 'nshf');
    monthlyMean = squeeze(mean(nshfData, [1, 2, 5], 'omitnan'));
    annualMean = mean(reshape(monthlyMean, 12, []), 1, 'omitnan');
    startYear = 1990 - ncread(filePath, 'year', 1, 1) + 1;
    endYear = 2050 - ncread(filePath, 'year', 1, 1) + 1;
    annualNshfMeans(:, i) = annualMean(startYear:endYear);
end

correlations = zeros(numModels, 1);
figure;
hold on;
colors = lines(numModels);

for i = 1:numModels
    correlations(i) = corr(annualWindMeans(:, i), annualNshfMeans(:, i), 'Rows', 'complete');
    scatter(annualWindMeans(:, i), annualNshfMeans(:, i), 36, colors(i, :), 'filled');
end

legend({nshfFiles.name});
title('Scatter Plot of Annual Mean Wind Speed vs. NSHF for Each Model');
xlabel('Annual Mean Wind Speed (m/s)');
ylabel('Annual Mean NSHF (W/m^2)');
hold off;

%% grid nshf, wind correlation labrador sea

labrador_lon_range = [-65, -50]; 
labrador_lat_range = [52, 62];   


labrador_lon_indices = find(lonRange >= labrador_lon_range(1) & lonRange <= labrador_lon_range(2));
labrador_lat_indices = find(latRange >= labrador_lat_range(1) & latRange <= labrador_lat_range(2));


labrador_dWind = squeeze(dWind(labrador_lon_indices, labrador_lat_indices, :));
labrador_dnshf = squeeze(dnshf(labrador_lon_indices, labrador_lat_indices, :));

labrador_correlation_matrix = zeros(length(labrador_lat_indices), length(labrador_lon_indices));

for lat = 1:length(labrador_lat_indices)
    for lon = 1:length(labrador_lon_indices)
        wind_changes = squeeze(labrador_dWind(lon, lat, :));
        nshf_changes = squeeze(labrador_dnshf(lon, lat, :));

        validIdx = ~isnan(wind_changes) & ~isnan(nshf_changes);
        if any(validIdx)
            labrador_correlation_matrix(lat, lon) = corr(wind_changes(validIdx), nshf_changes(validIdx));
        else
            labrador_correlation_matrix(lat, lon) = NaN;
        end
    end
end


labrador_lon_range_plot = lonRange(labrador_lon_indices);
labrador_lat_range_plot = latRange(labrador_lat_indices);

figure;
imagesc(labrador_lon_range_plot, labrador_lat_range_plot, labrador_correlation_matrix, 'AlphaData', ~isnan(labrador_correlation_matrix));
colorbar;
title('Correlation between Wind Speed Change and NSHF Change in Labrador Sea');
xlabel('Longitude');
ylabel('Latitude');
axis xy

%% all grid

entire_correlation_matrix = zeros(length(latRange), length(lonRange));

for lat = 1:length(latRange)
    for lon = 1:length(lonRange)
        wind_changes = squeeze(dWind(lon, end-lat+1, :)); 
        nshf_changes = squeeze(dnshf(lon, end-lat+1, :)); 

        validIdx = ~isnan(wind_changes) & ~isnan(nshf_changes);
        if any(validIdx)
            entire_correlation_matrix(lat, lon) = corr(wind_changes(validIdx), nshf_changes(validIdx));
        else
            entire_correlation_matrix(lat, lon) = NaN;
        end
    end
end

figure;
imagesc(lonRange, latRange, entire_correlation_matrix, 'AlphaData', ~isnan(entire_correlation_matrix));
colorbar;
title('Correlation between Wind Speed Change and NSHF Change (Entire Region)');
xlabel('Longitude');
ylabel('Latitude');

%%
annualWindMeans = zeros(length(years), numModels);
annualNshfMeans = zeros(length(years), numModels);

for i = 1:numModels
    filePath = fullfile(windFiles(i).folder, windFiles(i).name);
    windData = ncread(filePath, 'sfcWind');
    monthlyMeanWind = squeeze(mean(windData, [1, 2, 5], 'omitnan'));
    annualMeanWind = mean(reshape(monthlyMeanWind, 12, []), 1, 'omitnan');
    startYear = 1990 - ncread(filePath, 'year', 1, 1) + 1;
    endYear = 2050 - ncread(filePath, 'year', 1, 1) + 1;
    annualWindMeans(:, i) = annualMeanWind(startYear:endYear);
    
    filePath = fullfile(nshfFiles(i).folder, nshfFiles(i).name);
    nshfData = ncread(filePath, 'nshf');
    monthlyMeanNshf = squeeze(mean(nshfData, [1, 2, 5], 'omitnan'));
    annualMeanNshf = mean(reshape(monthlyMeanNshf, 12, []), 1, 'omitnan');
    annualNshfMeans(:, i) = annualMeanNshf(startYear:endYear);
end

correlations = zeros(numModels, 1);
for i = 1:numModels
    correlations(i) = corr(annualWindMeans(:, i), annualNshfMeans(:, i), 'Rows', 'complete');
end

figure;
hold on;
colors = lines(numModels);

for i = 1:numModels
    % Scatter plot
    scatter(annualWindMeans(:, i), annualNshfMeans(:, i), 36, colors(i, :), 'filled');
    
    p = polyfit(annualWindMeans(:, i), annualNshfMeans(:, i), 1);
    x_fit = linspace(min(annualWindMeans(:, i)), max(annualWindMeans(:, i)), 100);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, 'Color', colors(i, :), 'LineWidth', 2);
    
    text(mean(x_fit), mean(y_fit), sprintf('y = %.2fx + %.2f', p(1), p(2)), 'Color', colors(i, :));
    text(max(x_fit), min(y_fit), sprintf('Correlation = %.2f', correlations(i)), 'Color', colors(i, :), 'HorizontalAlignment', 'right');
end

hold off;
legend({nshfFiles.name}, 'Interpreter', 'none', 'Location', 'best');
title('Scatter Plot of Annual Mean Wind Speed vs. NSHF for Each Model with Regression Lines');
xlabel('Annual Mean Wind Speed (m/s)');
ylabel('Annual Mean NSHF (W/m^2)');
grid on;

fprintf('Correlations between annual mean wind speeds and NSHF for each model:\n');
disp(correlations);
%%

significant_models = find(abs(correlations) >= 0.6);

fprintf('Significant models (correlation >= 0.6):\n');
disp(significant_models);

figure;
hold on;
for i = 1:length(significant_models)
    idx = significant_models(i);
    scatter(annualWindMeans(:, idx), annualNshfMeans(:, idx), 36, colors(idx, :), 'filled');
end

for i = 1:length(significant_models)
    idx = significant_models(i);
    p = polyfit(annualWindMeans(:, idx), annualNshfMeans(:, idx), 1);
    x_fit = linspace(min(annualWindMeans(:, idx)), max(annualWindMeans(:, idx)), 100);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, 'Color', colors(significant_models(i), :), 'LineWidth', 2);
    text(mean(x_fit), mean(y_fit), sprintf('y = %.2fx + %.2f', p(1), p(2)), 'Color', colors(significant_models(i), :));
    corr_coef = corr(annualWindMeans(:, idx), annualNshfMeans(:, idx), 'Rows', 'complete');
    text(max(x_fit), min(y_fit), sprintf('Correlation = %.2f', corr_coef), 'Color', colors(significant_models(i), :), 'HorizontalAlignment', 'right');
end

hold off;

legend({nshfFiles(significant_models).name}, 'Interpreter', 'none', 'Location', 'best');
title('Scatter Plot of Annual Mean Wind Speed vs. NSHF for Significant Models (Correlation >= 0.6) with Regression Lines and Correlation Coefficients');
xlabel('Annual Mean Wind Speed (m/s)');
ylabel('Annual Mean NSHF (W/m^2)');

%% dNshf v Wind, dAmoc v Wind

numModels = numel(amocFiles);  % Number of models
startYear = 1990;
endYear = 2050;
numYears = endYear - startYear + 1;
decadeLength = 10;

damoc = zeros(numModels, 1);
dnshf = zeros(121, 31, numModels);
dWind = zeros(121, 31, numModels);

for i = 1:numModels
    currentFile = amocFiles(i).name;
    amoc_data = load(fullfile(folderPath, currentFile));
    startDecadeMean = mean(amoc_data((startYear-1980):(startYear-1980+decadeLength-1)));
    endDecadeMean = mean(amoc_data((endYear-1980-decadeLength+1):(endYear-1980)));
    damoc(i) = endDecadeMean - startDecadeMean;  % Change in AMOC
end

for i = 1:numModels
    % NSHF
    currentFile = nshfFiles(i).name;
    nshfPath = fullfile(folderPath, currentFile);
    nshfData = ncread(nshfPath, 'nshf');

    startYearIndex = (startYear - 1980) * 12 + 1;
    endYearIndex = (endYear - 1980) * 12;

    startNshf = mean(nshfData(:, :, startYearIndex:startYearIndex + decadeLength*12 - 1), 3, 'omitnan');
    endNshf = mean(nshfData(:, :, endYearIndex - decadeLength*12 + 1:endYearIndex), 3, 'omitnan');
    dnshf(:, :, i) = endNshf - startNshf;

    windPath = fullfile(folderPath, windFiles(i).name);
    windData = ncread(windPath, 'sfcWind');
    startWind = mean(windData(:, :, startYearIndex:startYearIndex + decadeLength*12 - 1), 3, 'omitnan');
    endWind = mean(windData(:, :, endYearIndex - decadeLength*12 + 1:endYearIndex), 3, 'omitnan');
    dWind(:, :, i) = endWind - startWind;
end

labrador_lon_range = [-65, -50]; 
labrador_lat_range = [52, 62];   
lonRange = -180:2:60;
latRange = 30:2:90;

labrador_lon_indices = find(lonRange >= labrador_lon_range(1) & lonRange <= labrador_lon_range(2));
labrador_lat_indices = find(latRange >= labrador_lat_range(1) & latRange <= labrador_lat_range(2));

labrador_dWind = dWind(labrador_lon_indices, labrador_lat_indices, :);
labrador_dnshf = dnshf(labrador_lon_indices, labrador_lat_indices, :);

avg_dWind = squeeze(mean(labrador_dWind, [1, 2], 'omitnan'));
avg_dnshf = squeeze(mean(labrador_dnshf, [1, 2], 'omitnan'));

[r_nshf, p_nshf] = corr(avg_dWind, avg_dnshf, 'Rows', 'complete');
[r_amoc, p_amoc] = corr(avg_dWind, damoc, 'Rows', 'complete');

figure;
subplot(1, 2, 1);
scatter(avg_dWind, avg_dnshf, 'filled');
hold on;
p1 = polyfit(avg_dWind, avg_dnshf, 1);
x_fit = linspace(min(avg_dWind), max(avg_dWind), 100);
y_fit = polyval(p1, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
title(['Wind speed vs. NSHF (r = ' num2str(r_nshf) ')']);
xlabel('ΔWind speed (m/s)');
ylabel('ΔNSHF (W/m^2)');
grid on;

subplot(1, 2, 2);
scatter(avg_dWind, damoc, 'filled');
hold on;
p2 = polyfit(avg_dWind, damoc, 1);
x_fit = linspace(min(avg_dWind), max(avg_dWind), 100);
y_fit = polyval(p2, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
title(['Wind speed vs. AMOC (r = ' num2str(r_amoc) ')']);
xlabel('ΔWind speed (m/s)');
ylabel('ΔAMOC (Sv)');
grid on;

fprintf('The regression line for NSHF is y = %.2fx + %.2f\n', p1(1), p1(2));
fprintf('The regression line for AMOC is y = %.2fx + %.2f\n', p2(1), p2(2));

%% dWindspeed

numModels = numel(windFiles);
startYear = 1990;
endYear = 2050;

startIndex = (startYear - 1980) * 12 + 1;
endIndex = (endYear - 1980 + 1) * 12;

lon_wind = ncread(fullfile(folderPath, windFiles(1).name), 'lon');
lat_wind = ncread(fullfile(folderPath, windFiles(1).name), 'lat');
numLon = length(lon_wind);
numLat = length(lat_wind);
meanWind1990 = zeros(numLon, numLat);
meanWind2050 = zeros(numLon, numLat);

for i = 1:numModels
    windPath = fullfile(folderPath, windFiles(i).name);
    windData = ncread(windPath, 'sfcWind');
    
    startWind = mean(windData(:, :, startIndex:startIndex+11), 3, 'omitnan');
    endWind = mean(windData(:, :, endIndex-11:endIndex), 3, 'omitnan');
    
    meanWind1990 = meanWind1990 + startWind / numModels;
    meanWind2050 = meanWind2050 + endWind / numModels;
end

deltaWind = meanWind2050 - meanWind1990;

figure;
axesm('MapProjection', 'eqdcylin', 'Frame', 'on', 'Grid', 'on', ...
    'MapLatLimit', [-90 90], 'MapLonLimit', [-180 180], ...
    'ParallelLabel', 'on', 'MeridianLabel', 'on');
surfm(lat_wind, lon_wind, deltaWind');
caxis([-0.3 0.3]);
colorbar('Ticks', -0.3:0.1:0.3, 'TickLabels', {'-0.3', '-0.2', '-0.1', '0', '0.1', '0.2', '0.3'});
colormap(parula);
title('\DeltaWind speed (2050-1990)', 'FontSize', 14);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow([land.Lat], [land.Lon], 'Color', 'k');