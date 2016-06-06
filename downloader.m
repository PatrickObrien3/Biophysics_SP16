% DOWNLOADER Download 221 PDB files to a subfolder named 'pdbs', for future
% use in Homework 4.
% 
% The subfolder 'pdbs' will be created if it does not already exist.

chains = readtable('chains.csv');
[N, ~] = size(chains);

%%
[~,~,~] = mkdir('pdbs');
for i=1:N
    strucname = chains{i, 1}{1};
    chainname = chains{i, 2}{1};
    struc = getpdb(strucname);
    fname = sprintf('pdbs/%s.pdb', strucname);
    struc = getpdb(strucname, 'ToFile', fname);
    disp(sprintf('%4d: %s', N-i, strucname));
end