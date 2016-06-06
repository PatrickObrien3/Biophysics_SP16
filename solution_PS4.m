
% Biophysics PS4
clear
chains = readtable('chains.csv');
[N, ~] = size(chains);


% Set the format of figures
figure(1);
clf;
hold on; box on; axis equal;
set(gca, 'linewidth', 2, 'fontsize', 20);
xlabel('\chi_1', 'fontsize', 20);
ylabel('\chi_2', 'fontsize', 20);
xlim([0, 360]); ylim([0, 360]);

figure(2);
clf;
hold on; box on; axis equal;
set(gca, 'linewidth', 2, 'fontsize', 20);
xlabel('\phi', 'fontsize', 20);
ylabel('\psi', 'fontsize', 20);
xlim([-180, 180]); ylim([-180, 180]);
ramachandranLimits();

for i=[1:82, 84:N]
    
       
    % Read the pdb file
    strucname = ['pdbs/', chains{i, 1}{1}, '.pdb'];
    struc = pdbread(strucname);

    % Fetch the useful data
    x = [struc.Model.Atom.X];
    y = [struc.Model.Atom.Y];
    z = [struc.Model.Atom.Z];
    resseq = [struc.Model.Atom.resSeq];
    resname = {struc.Model.Atom.resName};
    atomname = {struc.Model.Atom.AtomName};
    chainname = {struc.Model.Atom.chainID};
    
    % For some chains, the ResSeq doesn't start from 0
    if resseq(1) < 1
        resseq = resseq+1-resseq(1);
    end
    
    % Obtain the number of residues and the number of atoms
    N_res = max(resseq);
    [~, N_atom] = size(atomname);
    fprintf('%d, %d\n', i, N_atom);
    
    % Indicate which residues we have read
    % '1' means we have read the residual, '0' means not 
    mark_read = zeros(N_res, 1);
    
    % Store all \chi1 and \chi2 of LEU
    chi1 = zeros(N_res, 1);
    chi2 = zeros(N_res, 1);
    
    % Store all \phi and \psi
    phi = zeros(N_res, 1);
    psi = zeros(N_res, 1);
    
    % Store the positions of N, CA, and C in last residue
    N_last = 0;
    CA_last = 0;
    C_last = 0;
    
    for ii = 1:N_atom
        
        % Find the first atom of each residue
        if mark_read(resseq(ii)) == 0 
            mark_read(resseq(ii)) = 1;
            
            %% Calculate the side-chain dihedral angles for LEU
            if strcmp(chainname{ii}, chains{i, 2}{1}) && strcmp(resname{ii},'LEU') && strcmp(atomname{ii},'N') && strcmp(atomname{ii+1},'CA') && strcmp(atomname{ii+4},'CB') && strcmp(atomname{ii+5},'CG') && strcmp(atomname{ii+6},'CD1')
                % Specify the the position of N, CA, CB, and CG
                xyz1 = [x(ii), x(ii+1), x(ii+4), x(ii+5); y(ii), y(ii+1), y(ii+4), y(ii+5); z(ii), z(ii+1), z(ii+4), z(ii+5)];
                % Specify the the position of N, CA, CB, and CG
                xyz2 = [x(ii+1), x(ii+4), x(ii+5), x(ii+6); y(ii+1), y(ii+4), y(ii+5), y(ii+6); z(ii+1), z(ii+4), z(ii+5), z(ii+6)];
            
                % Calculate the \chi1 and \chi2
                chi1(resseq(ii)) = dihedral(xyz1);
                chi2(resseq(ii)) = dihedral(xyz2);
            end
            
            %% Calculate the backbone dihedral angels
            if ii+2>N_atom || strcmp(resname{ii},'PRO') || strcmp(resname{ii},'GLY') || ~strcmp(atomname{ii},'N') || ~strcmp(atomname{ii+1},'CA') || ~strcmp(atomname{ii+2},'C')
                N_last = 0;
                CA_last = 0;
                C_last = 0;
            else
                if N_last ~= 0 && resseq(ii) == resseq(N_last)+1
                    % Specify the the position of C, N, CA, and C
                    xyz3 = [x(C_last), x(ii), x(ii+1), x(ii+2); y(C_last), y(ii), y(ii+1), y(ii+2); z(C_last), z(ii), z(ii+1), z(ii+2)];
                    % Specify the the position of N, CA, C, and N
                    xyz4 = [x(N_last), x(CA_last), x(C_last), x(ii); y(N_last), y(CA_last), y(C_last), y(ii); z(N_last), z(CA_last), z(C_last), z(ii)];
                    
                    % Calculate the \phi and \psi
                    phi(resseq(ii)) = dihedral(xyz3);
                    psi(resseq(ii)) = dihedral(xyz4);

                end
                N_last = ii;
                CA_last =ii+1;
                C_last = ii+2;
                
            end
        end
        
    end
    
    % Plot the \chi2 vs \chi1 for LEU
    chi1(chi1==0) = [];
    chi1 = mod(chi1*180/pi, 360);
    chi2(chi2==0) = [];
    chi2 = mod(chi2*180/pi, 360);
    figure(1);
    plot(chi1, chi2, 'k.', 'markersize', 5);
    
    % Get the Ramachandran plot
    phi(phi==0) = [];
    phi = phi*180/pi;
    psi(psi==0) = [];
    psi = psi*180/pi;
    figure(2);
    plot(phi(1:end-1), psi(2:end), 'k.', 'markersize', 5);
    
    drawnow;    
end
