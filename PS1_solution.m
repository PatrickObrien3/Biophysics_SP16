% solution for PS1 in Biological Physics

clear;
clc;

%% A) A 1-Dimensional Random Walk
n_steps = 8;
n_tracks = 1000;
steps = 2*round(rand(n_steps,n_tracks))-1;
locations = cumsum(steps, 1);

figure(1); clf; hold on; box on;
set(gca, 'fontsize', 20);
set(gca, 'linewidth', 2);
plot(locations(:,1:6), 'linewidth', 2);
xlabel('Number of Steps', 'fontsize', 20);
ylabel('Location', 'fontsize', 20);

%% B) Statistics
final_locs = locations(n_steps,:);
mean_loc = mean(final_locs);
vari_loc = var(final_locs);

% mean_theory = p1*x1+p2*x2 = 0, 
% var_theory = N = 8;
mean_theory = 0;
var_theory = n_steps;

figure(2); clf; hold on; box on;
set(gca, 'fontsize', 20);
set(gca, 'linewidth', 2);

% Plot the histogram
bins = -n_steps-1:2:n_steps+1;
centers = -n_steps:2:n_steps;
histogram(final_locs,bins,'normalization','pdf');

% Calculate the binomial distribution
nl = 8:-1:0;
nr = 0:1:8;
bin_width = 2;
binom = factorial(n_steps) ./ (factorial(nl) .* factorial(nr)) .* 0.5^n_steps/bin_width;
plot(centers, binom,'ko-', 'linewidth', 2);

% Calculate the normal distribution
x_norm = -8:0.01:8;
y_norm = (1/sqrt(2*pi*var_theory)) .* exp(-(x_norm-mean_theory).^2 ./ (2 * var_theory));
plot(x_norm, y_norm, 'r-', 'linewidth', 2);
xlabel('Final Location', 'fontsize', 20);
ylabel('PDF', 'fontsize', 20);
title('1000 8-Step R.W. Histogram vs. Binomial/Gaussian', 'fontsize', 20);
legend('Histogram of Simulation','Binomial Distribution','Gaussian Distribution');
hold off;

%% Repeat aboove steps for n_steps = 1000 and n_tracks = 20000
n_steps = 1000; n_tracks = 20000;
mean_theory = 0;
var_theory = n_steps;
steps = 2*round(rand(n_steps,n_tracks))-1;
locations = cumsum(steps);

figure(3); 
clf; hold on; box on;
set(gca, 'fontsize', 20);
set(gca, 'linewidth', 2);
plot(locations(:,1:6), 'linewidth', 2);
xlabel('Number of Steps', 'fontsize', 20);
ylabel('Location', 'fontsize', 20);

figure(4); 
clf; hold on; box on;
set(gca, 'fontsize', 20);
set(gca, 'linewidth', 2);
final_locs = locations(n_steps,:);
bins = -100-5:10:100+5;
histogram(final_locs,bins,'normalization','pdf');

x_norm = -100:0.1:100;
y_norm = (1/sqrt(2*pi*var_theory)) .* exp(-(x_norm-mean_theory).^2 ./ (2*var_theory));
plot(x_norm, y_norm, 'r-', 'linewidth', 2);
xlabel('Final Location', 'fontsize', 20);
ylabel('PDF', 'fontsize', 20);
title('20000 1000-Step R.W. Histogram vs. Gaussian', 'fontsize', 20);
legend('Histogram of Simulation','Gaussian Distribution');


