% Biophysics PS2
% Please copy this file to the folder containing 'PS2data'

clear;
clc;

%% Load Data
load('PS2data');

dts = [1,100, 10e3]; % units of us
arrs = {locs1us, locs100us, locs10ms}; % units of um

%% Plot Data
figure(1);
clf;
for i=1:3
    arr = arrs{i};
    plot(arr(:,:,1), arr(:,:,2), 'k.', 'MarkerSize', .1);
    hold on;
end

cols = 'rbg';
for i=3:-1:1
    arr = arrs{i};
    plot(arr(:,1,1), arr(:,1,2), [cols(i) '-'], 'LineWidth', .1);
    hold on;
end

axis equal;
axis off;
saveas(gcf, 'cell-trajectories.png');

%% Calculate MSDs
figure(1); clf; hold on; box on;
set(gca, 'fontsize', 20, 'linewidth', 2, 'xscale', 'log', 'yscale', 'log');
xlabel('\Deltat (\mus)', 'fontsize', 20);
ylabel('MSD (\mum^2)', 'fontsize', 20);

MSDs3 = {};
ns3 = {};
for i=1:3
    arr = arrs{i};
    [tlen, Np, Ndim] = size(arr);
    ixs = unique(round(logspace(0, log10(tlen-1), 200)));
    MSDs = zeros(size(ixs));
    for j = 1:length(ixs)
        ix = ixs(j);
        dr = arr(1:end-ix,:,:) - arr(1+ix:end, :, :);
        drmag = dr(:,:,1).^2+dr(:,:,2).^2;
        %drmag = sum(dr.^2, 3);
        MSDs(j) = mean2(drmag);
    end
    
    plot(ixs * dts(i), MSDs, 'linewidth', 2);
    
    ns3{i} = ixs;
    MSDs3{i} = MSDs;
end

hgexport(gcf, 'three-MSDs.eps');

%% Calculate The Diffusion Coefficient and the Scaling Exponent
ix1 = ns3{1} * dts(1) < min(ns3{2}) * dts(2);
ix2 = ns3{2} * dts(2) < min(ns3{3}) * dts(3);

alldts = [ns3{1}(ix1)*dts(1), ns3{2}(ix2)*dts(2), ns3{3}*dts(3)];
allMSDs = [MSDs3{1}(ix1), MSDs3{2}(ix2), MSDs3{3}];

figure(3); clf; hold on; box on;
set(gca, 'fontsize', 20, 'linewidth', 2, 'xscale', 'log', 'yscale', 'log');
xlabel('\Deltat (\mus)', 'fontsize', 20);
ylabel('MSD (\mum^2)', 'fontsize', 20);

plot(alldts, allMSDs, 'linewidth', 2);

% Diffussion
ix = alldts < 10;
D = mean(allMSDs(ix) ./ alldts(ix)) / 4; % D=7.7337e-06 um^2/us
plot(alldts, 4*D*alldts, 'k:', 'linewidth', 2);

% Sub-diffusive
ix = (50 < alldts) & (alldts < 1000);
Y = log(allMSDs(ix));
X = log(alldts(ix));
fits = polyfit(X, Y, 1);
expt = fits(1);
C = exp(fits(2)); % C=7.4377e-05, expt=0.7008
plot(alldts, C .* (alldts.^expt), 'k--', 'linewidth', 2');

legend('MSD', 'Diffusion', 'Sub-diffusive', 'location', 'southeast');

hgexport(gcf, 'MSD-fit.eps');


