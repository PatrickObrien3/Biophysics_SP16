%% Biological Physics Problem Set 0
% Patrick O'Brien
% Due January 28th, 2016

%% Comparing Histograms

% random 1000 points with mean of 1 
ar1 = exprnd(1, [1000,1]);

% random 1000 points with mean of 3 
ar2 = exprnd(3, [1000,1]);

% random 2000 points with mean of 3 
ar3 = exprnd(3, [2000,1]);


histogram(ar1, 'normalization', 'pdf')
hold all
histogram(ar2, 'normalization', 'pdf')
hold all
histogram(ar3, 'normalization', 'pdf')
legend('Array 1: N=1000 u=1', 'Array 2: N=1000 u=3', 'Array 3: N=2000 u=3');

% Here I use the pdf normalization. 


%% Calculating a histogram
%clear all
rand_nums = exprnd(1, [1000,1]);
bins = linspace(0,5,21);
rand_counts = histcounts(rand_nums, bins);
figure
histogram(rand_nums, bins, 'normalization', 'pdf')
hold all
bin_centers = zeros(1, length(bins)-1);
for i = 1:length(bin_centers)
    bin_centers(i) = (bins(i) + bins(i+1))./(2);  
end

plot(bin_centers, exp(-bin_centers), 'linewidth', 2);
%title('Theoretical Distribution')


rand_counts2 = (rand_counts)./(length(rand_nums).*(diff(bins)));
plot(bin_centers, rand_counts2, 'linewidth',2 )
legend('histogram', 'theoretical distribution', 'rand counts')






