%% Biophysics Problem Set 05
%% Author: Patrick O'Brien
%% Date: April 21st, 2016

clear all
close all

%% Parameters

dt = 0.1;
r0 = 1;
rc = 0.5;
rho = 1.5;
L=6;
beta = [0,100];
v0 = 0.5;
etavals = [0.50:0.02:0.7];
Nsteps = 200;
nbins = 60;
edges = linspace(0,L/2,nbins);
delta = (edges(2)-edges(1));
N = rho * L * L;

%% Molecular Dynamics

% Basic setup
rs = rand(N, 2) * L;
vs = randn(N, 2);
vnorm = sqrt(sum(vs'.^2))';
vs = vs .* v0 ./ [vnorm, vnorm];
phis = zeros(length(etavals), Nsteps);

for k=1:length(beta)
    betaval = beta(k);
    for i=1:length(etavals)
        eta = etavals(i);
        for j = 1:Nsteps
            [vs] = vicsekvelocityJPO(v0, r0, eta, L, rs, vs, betaval, rc);
            rs = rs + vs * dt;
            phis(i,j) = sqrt(sum(mean(vs).^2)) ./ v0;
        end
    end
    phiavg(:,k) = mean(phis,2);
end

%% Plot of order parameter vs. eta 
% Below is the plot of the order parameter vs eta for non-repulsive and
% repulsive models. The repulsive graph has beta factor = 100. The value of
% eta_t, the value for which the order parameter is = 0.5, is shifted to a
% smaller value of eta when repulsive interactions are considered. 

figure, hold on
axis square
plot(etavals,phiavg,'LineWidth',1.5)
set(gca,'box','on');
legend('\beta = 0','\beta = 100');
xlabel('{\boldmath$\eta$}','Interpreter','latex','FontSize',13);
ylabel('{\boldmath$\left<\phi\right>$}','Interpreter','latex','FontSize',13);

%% Part B: The radial distribution function

etavals = [0.2,0.5,0.7];

for k=1:length(beta)
    betaval = beta(k);
    for i=1:length(etavals)
        num_sum = NaN(1,nbins-1);
        eta = etavals(i);
        for j = 1:Nsteps
            [vs,rijdists] = vicsekvelocityJPO(v0, r0, eta, L, rs, vs, betaval, rc);
            rs = rs + vs * dt;

            rijdists(rijdists == 0) = nan;
            rijdists(rijdists> L/2) = nan;

            [ns] = histcounts(rijdists,edges);
            num_sum = nansum([num_sum;ns],1);
        end
        ns_all(:,i,k) = num_sum;
    end
end

Npart = rho*L^2;
areas = 2*pi*edges(2:end)*delta;
norm = (Npart*Nsteps).*rho.*repmat(areas',[1 length(etavals) length(beta)]);
gr = ns_all./norm;
r_rc = edges(2:end)'./rc;

%% Plot Radial Distribution Functios g(r)

hardsphere = csvread('gr0.3.csv');
hsx = hardsphere(:,1);
hsy = hardsphere(:,2);

figure
hold on 
axis square
plot(r_rc,gr(:,:,1))
plot(hsx,hsy);
set(gca,'box','on');
legend('\eta = 0.2','\eta = 0.5','\eta = 0.7','Hard Sphere');
xlabel('r/r_{c}')
ylabel('g(r)')
title('Radial Distribution Function: \beta = 0')

figure
hold on 
axis square
plot(r_rc,gr(:,:,2))
plot(hsx,hsy);
set(gca,'box','on');
legend('\eta = 0.2','\eta = 0.5','\eta = 0.7','Hard Sphere');
xlabel('r/r_{c}')
ylabel('g(r)')
title('Radial Distribution Function: \beta = 100')