clear;
close all;
clc;

imagePath = '/Users/chikamaduabuchi/Dropbox (MIT)/RED LAB/MarcoProject/Segmentation/Frames/LAr_Baseline/1200.tif';

image = imread(imagePath);
binaryMask = imbinarize(image);

labeledBubbles = bwlabel(binaryMask);
bubbleMeasurements = regionprops(labeledBubbles, 'Perimeter', 'Area');
totalBubbles = length(bubbleMeasurements); 
fprintf('Total number of bubbles: %d\n', totalBubbles);

conversionFactor = 12.6;
perimeters_um = [bubbleMeasurements.Perimeter] * conversionFactor;
areas_um2 = [bubbleMeasurements.Area] * conversionFactor^2;
radii_um = sqrt(areas_um2 / pi);
radiiFromPerimeter_um = perimeters_um / (2 * pi);

%% Calculating the Dry Area and Contact Line
[d_daf, d_cld] = calc_density_fraction(binaryMask, conversionFactor);
fprintf('Dry Area Fraction (DAF): %f\n', d_daf);
fprintf('Contact Line Density (CLD): %f\n', d_cld);

%% Plotting Perimeter, Radii, and Area PDF/Count Distributions
plotDistributions(perimeters_um, areas_um2, radii_um);

%% Plotting the Radii from Perimeter and Area
plotRadiiComparison(radiiFromArea_um, radiiFromPerimeter_um);
plotPerimeterVsAreaRadius(radiiFromPerimeter_um, radiiFromArea_um);

%% Printing the Top 10 Occuring Radii
numBins = 15;
binEdges = linspace(min(radii_um), max(radii_um), numBins+1);
[~, ~, binIdx] = histcounts(radii_um, binEdges);
binCounts = accumarray(binIdx', 1, [numBins 1], @sum, 0);
[sortedCounts, sortIdx] = sort(binCounts, 'descend');
top10BinCounts = sortedCounts(1:10);
top10BinEdges = binEdges(sortIdx(1:10));

binStarts = top10BinEdges;
binEnds = binEdges(sortIdx(1:10) + 1);

top10BinsTable = table(binStarts', binEnds', top10BinCounts, 'VariableNames', {'BinStart_um', 'BinEnd_um', 'Frequency'});

disp('Top 10 most occurring radii bins (um) and their frequencies:');
disp(top10BinsTable);
filename = 'top10BinsTable.csv';
writetable(top10BinsTable, filename);
disp(['Top 10 bins table written to ' filename]);


%% Local Functions
function plotDistributions(perimeters, areas, radii)
    plotAttributeDistributions(perimeters, 'Perimeter (µm)');
    plotAttributeDistributions(areas, 'Area (µm^2)');
    plotAttributeDistributions(radii, 'Radius (µm)');
end

function plotAttributeDistributions(data, dataLabel)
    numBins = 15;
    binsLog = createLogBins(data, numBins);
    figure;
    subplot(2, 2, 1);
    plotPDF(data, 'Linear', [], 'linear', 'pdf', dataLabel);
    subplot(2, 2, 2);
    plotPDF(data, 'Log', binsLog, 'log', 'pdf', dataLabel);
    subplot(2, 2, 3);
    plotPDF(data, 'Linear', [], 'linear', 'count', dataLabel);
    subplot(2, 2, 4);
    plotPDF(data, 'Log', binsLog, 'log', 'count', dataLabel);
end

function binsLog = createLogBins(data, numBins)
    minData = min(data(data > 0));
    maxData = max(data);
    binsLog = logspace(log10(minData), log10(maxData), numBins);
end

function plotPDF(data, plotTitle, bins, axisType, normalization, dataLabel)
    if isempty(bins)
        histogram(data, 'Normalization', normalization);
    else
        histogram(data, 'BinEdges', bins, 'Normalization', normalization);
    end
    set(gca, 'XScale', axisType);
    if strcmp(axisType, 'log')
        minData = min(data(data > 0));
        maxData = max(data);
        ticks = logspace(floor(log10(minData)), ceil(log10(maxData)), ceil(log10(maxData))-floor(log10(minData))+1);
        set(gca, 'XTick', ticks);
        tickLabels = arrayfun(@(v) sprintf('10^{%d}', round(log10(v))), ticks, 'UniformOutput', false);
        set(gca, 'XTickLabel', tickLabels);
    end
    xlabel(dataLabel);
    if strcmp(normalization, 'pdf')
        if contains(dataLabel, 'Area')
            ylabel('Probability Density (1/\mum^2)');
        elseif contains(dataLabel, 'Perimeter')
            ylabel('Probability Density (1/\mum)');
        else
            ylabel('Probability Density (1/\mum)');
        end
    else
        ylabel('Count');
    end
    title(plotTitle);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'FontSize', 12);
    set(gca, 'FontWeight', 'bold');
    grid on;
    set(gca, 'GridLineStyle', '-');
    set(gca, 'MinorGridLineStyle', ':');
    set(gca, 'GridAlpha', 0.2);
    set(gca, 'MinorGridAlpha', 0.2);
end

function plotPerimeterVsAreaRadius(radiiFromPerimeter, radiiFromArea)
    figure;
    scatter(radiiFromPerimeter, radiiFromArea, 'filled');
    xlabel('Radius from Perimeter (µm)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Radius from Area (µm)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Scatter Plot of Radii: Perimeter vs. Area', 'FontSize', 14, 'FontWeight', 'bold');
    set(gca, 'LineWidth', 2, 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'GridAlpha', 0.5, 'MinorGridAlpha', 0.2, 'Box', 'on'); % Adjust alpha for grid lines
    set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on');
    grid on;
end

function plotRadiiComparison(radiiFromArea, radiiFromPerimeter)
    figure;
    hold on;
    plot(radiiFromArea, 'r', 'LineWidth', 2);
    plot(radiiFromPerimeter, 'b--', 'LineWidth', 2);
    legend('Radius from Area', 'Radius from Perimeter', 'Location', 'best');
    xlabel('Bubble Index', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Radius (µm)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Comparison of Radii Calculated from Area and Perimeter', 'FontSize', 14, 'FontWeight', 'bold');
    set(gca, 'LineWidth', 2, 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'GridAlpha', 0.5, 'MinorGridAlpha', 0.2, 'Box', 'on'); % Adjust alpha for grid lines
    set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on');
    grid on;
    hold off;
end

function [d_daf, d_cld] = calc_density_fraction(binary_mask, conversionFactor)
    pixelArea_um2 = conversionFactor^2; 
    total_pixels = numel(binary_mask);
    wet_pixels = sum(binary_mask == 0, 'all');
    d_daf = 1 - (wet_pixels / total_pixels);

    inverted_binary_mask = ~binary_mask;
    dist = bwdist(inverted_binary_mask, 'euclidean');
    dist(dist > 1) = 0;
    contact_line_length_pixels = sum(dist(:));
    contact_line_length_um = contact_line_length_pixels * conversionFactor;
    total_area_um2 = total_pixels * pixelArea_um2;
    d_cld = contact_line_length_um / total_area_um2;
end