% solution for PS0 in Biological Physics

clear;
clc;

%% Comparing Histograms

rand1 = exprnd(1,1,1000);
rand2 = exprnd(3,1,1000);
rand3 = exprnd(3,1,2000);

figure(1);
clf; hold on; box on;
set(gca, 'fontsize', 20, 'linewidth', 2);
title('Histograms from Exponential Distributions');
xlabel('x', 'fontsize', 20);
ylabel('PDF(x)', 'fontsize', 20);

histogram(rand1, 20, 'normalization', 'pdf');
histogram(rand2, 20, 'normalization', 'pdf');
histogram(rand3,20, 'normalization', 'pdf');
% We use 'pdf' to normalize, because the 'pdf' would gurantee the sum of
% the area of all bars is 1, which makes the histograms from different distribution comparable 

legend('mean=1, 1000points', 'mean=3, 1000points', 'mean=3, 2000points');

%% Calculating a histogram

rand_nums = exprnd(1,1,1000);
bins = linspace(0,5,21);
rand_counts = histcounts(rand_nums,bins);

figure(2);
clf; hold on; box on;
set(gca, 'fontsize', 20, 'linewidth', 2);
title('Histogram, Normalized Curve, and Theretical Curve');
xlabel('x', 'fontsize', 20);
ylabel('PDF(x)', 'fontsize', 20);

% plot the histogram
histogram(rand_nums,bins,'normalization','pdf');

% plot the theoretical curve
bin_centers = 0.5*(bins(1:end-1)+bins(2:end));
plot(bin_centers, exp(-bin_centers), '--', 'linewidth', 2);

% plot the normalized curve
rand_counts_norm = rand_counts./(sum(rand_counts).*diff(bins));
plot(bin_centers,rand_counts_norm,  'linewidth', 2);

legend('Histogram', 'Theoretical Curve', 'Normalized Curve');


