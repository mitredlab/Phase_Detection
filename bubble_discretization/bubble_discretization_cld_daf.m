clc; close all; clear;

L = 0.001;
numSimulations = 100;
N_values = 20:10:250;
R_values = linspace(10E-6, 100E-6, 10);
numN = length(N_values);
numR = length(R_values);

ME_Area = zeros(numSimulations, numN, numR);
PRE_Area = zeros(numSimulations, numN, numR);
ME_Perimeter = zeros(numSimulations, numN, numR);
PRE_Perimeter = zeros(numSimulations, numN, numR);
ME_DAF = zeros(numSimulations, numN, numR);
PRE_DAF = zeros(numSimulations, numN, numR);
ME_CLD = zeros(numSimulations, numN, numR);
PRE_CLD = zeros(numSimulations, numN, numR);

rng(0);
totalIterations = numSimulations * numN * numR;
currentIteration = 0;

tic;

for N_idx = 1:numN
    for R_idx = 1:numR
        N = N_values(N_idx);
        R = R_values(R_idx);
        [Xc, Yc] = meshgrid(linspace(L/N/2, L-L/N/2, N));
        A_theoretical = pi * R^2;
        P_theoretical = 2 * pi * R;
        A_total = L^2;
        DAF_Theoretical = A_theoretical / A_total;
        CLD_Theoretical = P_theoretical / A_total;

        for k = 1:numSimulations
            currentIteration = currentIteration + 1;
            Xb = R + (L - 2*R) * rand();
            Yb = R + (L - 2*R) * rand();
            Dis = sqrt((Xc - Xb).^2 + (Yc - Yb).^2);
            bub = (Dis < R);

            areaSimulated = sum(bub, 'all') * (L/N)^2;
            perimeterSimulated = sum(bwperim(bub), 'all') * (L/N);

            [d_daf, d_cld] = calc_density_fraction(bub);

            ME_Area(k, N_idx, R_idx) = A_theoretical - areaSimulated;
            PRE_Area(k, N_idx, R_idx) = 100 * (A_theoretical - areaSimulated) / A_theoretical;
            ME_Perimeter(k, N_idx, R_idx) = P_theoretical - perimeterSimulated;
            PRE_Perimeter(k, N_idx, R_idx) = 100 * (P_theoretical - perimeterSimulated) / P_theoretical;
            ME_DAF(k, N_idx, R_idx) = DAF_Theoretical - d_daf;
            PRE_DAF(k, N_idx, R_idx) = 100 * (DAF_Theoretical - d_daf) / DAF_Theoretical;
            ME_CLD(k, N_idx, R_idx) = CLD_Theoretical - d_cld;
            PRE_CLD(k, N_idx, R_idx) = 100 * (CLD_Theoretical - d_cld) / CLD_Theoretical;

            if mod(currentIteration, 100) == 0
                fprintf('Progress: %.2f%%\n', (currentIteration / totalIterations) * 100);
            end
        end
    end
end

totalTime = toc;
fprintf('Total simulation time: %.2f seconds.\n', totalTime);

save('bubble_discretization_results_cld_daf.mat', 'ME_Area', 'PRE_Area', 'ME_Perimeter', 'PRE_Perimeter', 'ME_DAF', 'PRE_DAF', 'ME_CLD', 'PRE_CLD');


%% Visuals
clc; close all;

% Load the results
load('simulation_errors.mat', 'ME_Area', 'PRE_Area', 'ME_Perimeter', 'PRE_Perimeter', 'ME_DAF', 'PRE_DAF', 'ME_CLD', 'PRE_CLD');

% Set figure properties
figPosition = [100, 100, 1024, 768];
fontSize = 20;
lineWidth = 2;
numBins = 20;

% Histograms for Area and Perimeter
fig1 = figure;
set(fig1, 'Position', figPosition);

subplot(2, 2, 1);
histogram(ME_Area(:), numBins);
xlabel('Mean Error of Area (\mum^2)', 'FontSize', fontSize);
ylabel('Frequency', 'FontSize', fontSize);
title('Mean Error for Area', 'FontSize', fontSize);
set(gca, 'LineWidth', lineWidth, 'FontSize', 16, 'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');

subplot(2, 2, 2);
histogram(ME_Perimeter(:), numBins);
xlabel('Mean Error of Perimeter (\mum)', 'FontSize', fontSize);
ylabel('Frequency', 'FontSize', fontSize);
title('Mean Error for Perimeter', 'FontSize', fontSize);
set(gca, 'LineWidth', lineWidth, 'FontSize', 16, 'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');

subplot(2, 2, 3);
histogram(PRE_Area(:), numBins);
xlabel('PRE of Area (%)', 'FontSize', fontSize);
ylabel('Frequency', 'FontSize', fontSize);
title('PRE for Area', 'FontSize', fontSize);
set(gca, 'LineWidth', lineWidth, 'FontSize', 16, 'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');

subplot(2, 2, 4);
histogram(PRE_Perimeter(:), numBins);
xlabel('PRE of Perimeter (%)', 'FontSize', fontSize);
ylabel('Frequency', 'FontSize', fontSize);
title('PRE for Perimeter', 'FontSize', fontSize);
set(gca, 'LineWidth', lineWidth, 'FontSize', 16, 'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');

% Histograms for DAF and CLD
fig2 = figure;
set(fig2, 'Position', figPosition);

subplot(2, 2, 1);
histogram(ME_DAF(:), numBins);
xlabel('Mean Error of DAF', 'FontSize', fontSize);
ylabel('Frequency', 'FontSize', fontSize);
title('Mean Error for DAF', 'FontSize', fontSize);
set(gca, 'LineWidth', lineWidth, 'FontSize', 16, 'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');

subplot(2, 2, 2);
histogram(ME_CLD(:), numBins);
xlabel('Mean Error of CLD', 'FontSize', fontSize);
ylabel('Frequency', 'FontSize', fontSize);
title('Mean Error for CLD', 'FontSize', fontSize);
set(gca, 'LineWidth', lineWidth, 'FontSize', 16, 'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');

subplot(2, 2, 3);
histogram(PRE_DAF(:), numBins);
xlabel('PRE of DAF (%)', 'FontSize', fontSize);
ylabel('Frequency', 'FontSize', fontSize);
title('PRE for DAF', 'FontSize', fontSize);
set(gca, 'LineWidth', lineWidth, 'FontSize', 16, 'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');

subplot(2, 2, 4);
histogram(PRE_CLD(:), numBins);
xlabel('PRE of CLD (%)', 'FontSize', fontSize);
ylabel('Frequency', 'FontSize', fontSize);
title('PRE for CLD', 'FontSize', fontSize);
set(gca, 'LineWidth', lineWidth, 'FontSize', 16, 'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');
%%
% Local Functions
function [d_daf, d_cld] = calc_density_fraction(binary_mask)
    total_pixels = numel(binary_mask);
    wet_pixels = sum(binary_mask == 0, 'all');
    d_daf = (1 - (wet_pixels / total_pixels));

    inverted_binary_mask = ~binary_mask;
    dist = bwdist(inverted_binary_mask);
    dist(dist > 1) = 0;
    contact_line_length = sum(dist, 'all');
    d_cld = contact_line_length / total_pixels;
end
