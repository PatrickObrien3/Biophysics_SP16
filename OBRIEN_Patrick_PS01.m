%% Biological Physics Problem Set 01
%% Patrick O'Brien
%% Due Thurs., Feb. 4, 2016

%% A 1-Dimensional Random Walk

% Produce a 8x1000 array of steps, where each is -1 or +1
clear
close all

steps = -1+(2).*(round(rand(8,1000)));

% Now make an array of locations using cumsum. Want 1000 tracks of 8 steps
% each, therefore will want to sum the columns. 
figure
locations = cumsum(steps);
plot(locations(:,1:6))
xlabel('Step Number');
ylabel('Position');
title('Position vs. Step Number for 1D Random Walk; N = 8')


%% Statistics

% Use indexing to get a 1x1000 array of final location of each track

floc = locations(8,:);

% calculated mean:
floc_u = mean(floc);
disp(floc_u,'calculated mean')
% Theoretical mean = 0. 

% calculated variance:
floc_var = var(floc);
disp(floc_var,'calculated variance')

% Theoretical variance = N = 8. 

% Histogram of P(x). 
figure
bins = [-9:2:9];
histogram(floc, bins, 'normalization', 'pdf')

% binomial distribution

hold all
bin_centers = zeros(1, length(bins)-1);
for i = 1:length(bin_centers)
    bin_centers(i) = (bins(i) + bins(i+1))./(2);  
end
% P(nl,nr) = (nl + nr)!/nl!.*nr! .* pr.*pl 
p = 0.5;
N = 8; 
nr = (0:8);
binomial = (factorial(N)./(factorial(N-nr).*factorial(nr))).*(0.5).^(N)./2;
plot(bin_centers, binomial, 'ko-')

% Normal Distribution

n = 8; 
xvals = linspace(-10, 10,200);
yvals = (1./( sqrt(2.*pi.*n) )).*exp( (-xvals.^(2)) ./(2.*n));
plot(xvals, yvals, 'r-');
xlabel('x')
ylabel('P(x)')
legend('Simulation Distribution', 'Binomial Distribution', 'Normal Distribution')
title('Distributions of Displacement Probabilities with N = 8')

%% A 1-Dimensional Random Walk With 1000 steps

% Produce a 1000,20000 array of steps, where each is -1 or +1
clear
close all

steps = -1+(2).*(round(rand(1000,20000)));
 
figure
locations = cumsum(steps,1);
plot(locations(:,1:6))
xlabel('Step Number');
ylabel('Position');
title('Position vs. Step Number for 1D Random Walk; N = 1000')


%% Statistics 

% Use indexing to get a 1x1000 array of final location of each track
floc = locations(end,:);

% Histogram of P(x). 
figure
bins = linspace(-120,120, 61);
histogram(floc, bins, 'normalization', 'pdf')
hold all

n = 1000; 
xvals = linspace(-100, 100,200);
yvals = (1./(sqrt(2.*pi.*n))).*exp((-xvals.^(2))./(2.*n));
plot(xvals, yvals, 'r-');
xlabel('x')
ylabel('P(x)')
legend('Simulation Distribution', 'Normal Distribution')
title('Distributions of Displacement Probabilities with N = 1000')











