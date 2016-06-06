%% Biological Physics Problem Set 03
%% Patrick O'Brien
%% Due Thurs., March. 3, 2016

clear
close all


%% A) Creating Chains

N = 200; % number of monomers
m = 500; % number of simulations per polymer
b = 1; % bond length
angles = 2.*pi.*rand((N-1),m); % random angles betweeen 0,2pi

dx = b.*cos(angles);
dy = b.*sin(angles);

x = cumsum(dx);
x = [zeros(1,500); x];

y = cumsum(dy);
y = [zeros(1,500); y];

%% plot first chain: 
 
figure(1)
plot(x(:,1), y(:,1))
axis equal
title('first chain of 200 monomers')
xlabel('x displacement')
ylabel('y displacement')

%% B) FJC Statistics

% Generate end vectors:
% xend = x(end,:); % Last row of each x and y, since starting from (0,0).
% yend = y(end,:); % Gets 1x500, expected, 500 different chains. 

endvec = [x(end,:); y(end,:)]; % produces a 2x500 matrix, so each column is a x,y end vector
% The average end-to-end vector is given here: 
endavg = mean(endvec,2);
disp(['Average end-to-end vector <R> = [', num2str(endavg(1)),', ', num2str(endavg(2)), ']'])
% Now R^2
R2 = mean(endvec.^2,2);
% text2 = ['Average end-to-end vector squared is: x = ', num2str(R2(1,1)), ' y = ', num2str(R2(2,1))];
% disp(text2)
% If I sum the two components they equal about 200, i.e. 187? 
% R2len = sqrt((R2(1,1).^2) + (R2(2,1).^2));
R2len = R2(1,1) + R2(2,1);
R2tlen = N.*b.^2;
disp(['Calculated <R^2>: ', num2str(R2len)]);
disp(['Theoretical <R^2>: ',num2str(R2tlen)]);

%% Vector of magnitude of endvec

Rmag = sqrt(endvec(1,:).^2 + endvec(2,:).^2);
[n,xout] = hist(Rmag, 18); % xout hold x coordinates (magnitudes) and the n
% holds the number of chains that have that value of xout

Pt = (2.*Rmag./N).*exp(-((Rmag.^2)./N)); % Theoretical distribution
binw = xout(2)-xout(1); %binwidth
area = cumsum(n.*binw); % creates vector with cumsum 
nnorm = n./area(1,end);

figure % Fig 2 stairs
stairs(xout, nnorm)
hold on
plot(Rmag,Pt, 'k.')
title('End-to-End Distance for the 2D FJC Model, N=200')
xlabel('|R|')
ylabel('P(|R|)')
legend('|R| simulated, 500 runs', '|R| theoretical')

%% C) Plotting a Single Chain

sepvecs = frc(pi./8, 200, 1); 
figure(3)
plot3(cumsum(sepvecs(1,:)),cumsum(sepvecs(2,:)),cumsum(sepvecs(3,:)),'.-')
axis equal
title('3D plot of 1 chain from FRC model with N = 200 monomers')
xlabel('x position')
ylabel('y position')
zlabel('z position')

%% D) Analysis of <ri rj>

sv = frc(pi./3, 200, 1000);
% make 200x1000 mag ri values 
rmag = sqrt(sv(1,:,:).^2 + sv(2,:,:).^2 + sv(3,:,:).^2);
disp(['The mean is: ', num2str(mean(mean(rmag))), ' and the stdev is: ', num2str(std(std(rmag)))])
disp('The mean should be = 1 because bond length = 1')
disp('The standard deviation is non-zero because of ')
disp('Matlab rounding errors in calculations')
rdot = zeros(1,199,1000);
for i = 1:199
    for j = 1:1000
        rdot(1,i,j) = dot(sv(:,i,j),sv(:,i+1,j));
    end
end
disp(['ri *r(i+1) mean is ', num2str(mean(mean(rdot))),' and the stdev is ',num2str(std(std(rdot)))])
% This matches theory, since 1.^2cos(pi/3) = 0.5. 
disp('The mean matches theory which is equal to')
disp('b^2*cos(theta), which here I calculate a ')
disp('value of 1^2cos(pi/3) = 0.5. The stdev   ')
disp('explanation is the same as above, from ')
disp('rounding erros in Matlab')
%% Plot ri dot rj vs i-j for theta = ...

theta_vals = [pi./100, pi./8, pi./3, pi./2, (2.*pi)./3];
lpcalc = zeros(1,length(theta_vals));

figure(4)
rstore = zeros(5,N);
for k = 1:length(theta_vals) 
    sv = frc(theta_vals(k), 200, 1000); 
    % array of d values; d = |i-j|
    dvals = zeros(1,N-1);
    for j = 1:N
        dvals(j) =abs(1-j);
    end
    for i = 1:length(dvals)
        M = sv(:,1:N-i,:).*sv(:,1+i:N,:);
        B = sum(M,1);
        rstore(k,i) = mean(mean(B));
    end
    plot(dvals,rstore(k,:),'.','LineWidth', 1)
    hold on
    rt = cos(theta_vals(k)).^dvals;
    plot(dvals, rt,':','LineWidth', 1)

end
title('<ri*rj> vs. |i-j|')  
xlabel('|i-j|')
ylabel('<ri*rj>')
legend('pi/100 sim', 'pi/100 theor', 'pi/8 sim', 'pi/8 theor', 'pi/3 sim', 'pi/3 theor', 'pi/2 sim', 'pi/2 theor', '2pi/3 sim', '2pi/3 theor') 
xlim([0 20])
ylim([-0.5 2])

%% Next section

disp('   calc lp    theoretical lp')
disp(['    ',num2str(find(rstore(1,:) > exp(-1),1,'last')),'       ',num2str(real(-1./log(cos(theta_vals(1)))))])      
disp(['    ',num2str(find(rstore(2,:) > exp(-1),1,'last')),'        ',num2str(real(-1./log(cos(theta_vals(2)))))])
disp(['    ',num2str(find(rstore(3,:) > exp(-1),1,'last')),'         ',num2str(real(-1./log(cos(theta_vals(3)))))])
disp(['    ',num2str(lpcalc(4)),'         ',num2str(real(-1./log(cos(theta_vals(4)))))])
disp(['    ',num2str(lpcalc(5)),'         ',num2str(real(-1./log(cos(theta_vals(5)))))])

% for each theta, calculate persistence length from data and also
% theoretical one, lp = -b./(ln(cos(theta)). 

% Explanation:
disp('From left to right, as the magnitude of i-j increases the monomers')
disp('get further apart. The <ri*rj> represents how correlated the motions')
disp('are. Small values of theta give more persistant polymers and a')
disp('higher persistence length since there is a further distance each ')
disp('polymer will extend before diverging from a straight path. That is ')
disp('why for pi/100, gives an almost infinity persistance length (and for')
disp('our simulation value it is as long as the number of monomers). There')
disp('is little change in monomer angle so the end-to-end stays in the')
disp('same direction. As theta increases, it is harder for the polymer to')
disp('be long. This is seen by the observation that as theta increases')
disp('the value of<ri*rj> goes to zero in less number of monomers away.')
%% E) sqrroot (R^2) vs. lp

lpvals = logspace(-2,4, 20); % 20 vals log spaced
R2 = zeros(1,20);
theta_vals = acos(exp(-1./lpvals));

for i = 1:length(lpvals)
    theta = theta_vals(i);
    sv = frc(theta,N, 1000);
    S = cumsum(sv,2);
    dist = sqrt(S(1,end,:).^2+ S(2,end,:).^2+ S(3,end,:).^2);
    R2(i) = sqrt(mean(dist.^2));
end

figure(5)
semilogx(lpvals, R2,'r.');
ylim([0 200])
hold all

lpvals2 = logspace(-2,4,200);
cosvals = exp(-1./lpvals2);

FRC_R2 = sqrt(N.*(( (1+ cosvals)./(1-cosvals)) - (2.*cosvals.*(1-(cosvals).^N))./(N.*(1-cosvals).^2)));
plot(lpvals2, FRC_R2, 'k-')

WLC_R2 = sqrt((2.*lpvals2.*N.*b)-2.*(lpvals2).^2.*(1-exp(-(N./lpvals2))));
plot(lpvals2, WLC_R2, 'b-')
legend('simulation', 'FRC', 'Worm-like Chain', 'location', 'northwest')
title('sqrt(<R^2>) for FRC, N =200')
xlabel('Persistence Length')
ylabel('sqrtR^2')

disp('The FRC is accurate for the entire simulation, whereas the worm-like')
disp('model is only accurate for persistence lengths greater than 1. The ')
disp('approximation used is that theta is close to zero, which would give')
disp('a high persistence length. Therefore, the WLC fails at small lp values')
disp('(or high theta values) which is seen in the graph.')







