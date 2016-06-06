%% Biophysics Problem Set 04
%% Patrick O'Brien
%% March 23, 2016

%% Part 1. Download Structures in chains
% This was done with downloader.m
%% Part 2. Calculate Leucine Dihedral angles and calculating phi and psi

clear 
close all

tic

chains = readtable('chains.csv');
[N,~] = size(chains);

X1 = [0];
X2 = [0];
phi = [0];
psi = [0];

for i = 1:N
    strucname = chains{i,1}{1};
    chainname = chains{i,2}{1};
    fname = sprintf('pdbs/%s.pdb', strucname);
    s = pdbread(fname);
    a = s.Model.Atom;
    M1 = struct2cell(a);
    M2 = reshape(M1,size(M1,1),size(M1,3))';
    M3 = M2;
    ina = strcmp(M2(:,3), '');
    inb = strcmp(M2(:,3), 'A');
    in0 = logical(ina+inb);
    M2 = M2(in0,:);
    
    sortchain = strcmp(M2(:,5),chainname);
    M2 = M2(sortchain,:);
    atomresnum = cell2mat(M2(:,6));
    totres = atomresnum(end); % total number of residues.
    
    in = strcmp(M2(:,4),'LEU');
    L3 = M2(in,:);                     % Get Leucines
    c6 = cell2mat(L3(:,6));
    aNin  = strcmp(L3(:,2),'N');
    aN = L3(aNin,:);                   % Get N
    aCAin = strcmp(L3(:,2),'CA');
    aCA = L3(aCAin,:);                 % Get CA
    aCBin = strcmp(L3(:,2),'CB');
    aCB = L3(aCBin,:);                 % Get CB
    aCGin = strcmp(L3(:,2),'CG');
    aCG = L3(aCGin,:);                 % Get CG
    aCDin = strcmp(L3(:,2),'CD1');
    aCD1 = L3(aCDin,:);                % Get CD
    numL = min([sum(aCAin) sum(aCGin) sum(aNin) sum(aCBin) sum(aCDin)]);
    % Use min to ensure gets done if multiple entries. 
    for k = 1:numL
        X1xyz = [aN(k,8:10); aCA(k,8:10); aCB(k,8:10); aCG(k,8:10)]';
        X2xyz = [aCA(k,8:10); aCB(k,8:10); aCG(k,8:10); aCD1(k,8:10)]';
        X1 = [X1 ; dihedral(cell2mat(X1xyz))];
        X2 = [X2 ; dihedral(cell2mat(X2xyz))];
    end

%     runningtotal = runningtotal + totres -2; 
    e = unique(cell2mat(M2(:,6)));
    for j = e(2):e(end-1) % This will exclude the first and last residue
        in = ismembc(atomresnum,j);
        res = M2(in,:);
        
        Nin = strcmp(res(:,2),'N');
        Cin = strcmp(res(:,2),'C');
        CAin = strcmp(res(:,2),'CA');
 
        inC = ismembc(atomresnum,j-1);
        resm1 = M2(inC,:);
        Cm1in = strcmp(resm1(:,2),'C');

        inN = ismembc(atomresnum,j+1);
        resp1 = M2(inN,:);
        Np1in = strcmp(resp1(:,2),'N');
        
        if isempty(res) == 0 && strcmp(res(1,4),'GLY') ~= 1 && strcmp(res(1,4),'PRO') ~= 1 ...
            && sum(Cm1in) ==1 && sum(Np1in) ==1 && sum(Nin) ==1 && sum(Cin) ==1 && sum(CAin) == 1
            
        Na = res(Nin,8:10);
        C = res(Cin,8:10);
        CA = res(CAin,8:10);
        Cm1 = resm1(Cm1in,8:10);
        Np1 = resp1(Np1in,8:10);
            
        x = dihedral(cell2mat([Cm1(1,:); Na(1,:); CA(1,:); C(1,:)]'));
        y = dihedral(cell2mat([Na(1,:); CA(1,:); C(1,:); Np1(1,:)]'));
        phi = [phi ; x];
        psi = [psi ; y];   
        end    
    end 
end

X1 = mod(X1*180/pi,360);
X2 = mod(X2*180/pi,360);

X1 = X1(2:end);
X2 = X2(2:end);

%% Plotting Leucine Dihedral Angles
figure
plot(X1,X2, 'k.')
xlabel('LEU Chi 1')
ylabel('LEU Chi 2')
xlim([0 360])
ylim([0 360])
title('Chi 1 vs. Chi 2 for Leucine')

%% Plotting the dihedral angles 
phi = phi(2:end);
psi = psi(2:end);

phi = rad2deg(phi);
psi = rad2deg(psi);

figure
plot(phi,psi, 'k.')
xlim([-180 180])
ylim([-180 180])
xlabel('Phi')
ylabel('Psi')
title('Ramachandran Plot for 221 pdb files')
ramachandranLimits()

toc
