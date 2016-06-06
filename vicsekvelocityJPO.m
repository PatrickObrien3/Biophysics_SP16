function [fs,rijdists] = vicsekvelocityJPO(v0, r0, eta, L, rs, vs, beta, rc)
% FORCES Finds the updated velocity for a vicsek model.
% From http://journals.aps.org/pre/abstract/10.1103/PhysRevE.77.046113#fulltext
% also at http://dx.doi.org/10.1103/PhysRevE.77.046113
%
% Returns a N x 2 array for the new velocities for each particle in the system.
    [N, ~] = size(rs);
    sum_vs = vs;
    Nis = ones(N, 1);
    fij = zeros(N, 2);
    fijtot = zeros(N, 2);
    rijdists = NaN(N,N);
    
    for i=1:N
       rijs = [rs(:, 1) - rs(i, 1), rs(:, 2) - rs(i, 2)];
       rijs = mod((rijs + L./2), L) - L./2;
       dists = sqrt(sum(rijs'.^2))';
       rijdists(:,i) = dists;
       in_range = dists < r0;
       rep_range = dists < rc;
       unitvec = normer(rijs);
       fij = ((1-rijs./rc).*unitvec).*[rep_range,rep_range];
       fijtot(i,:) = nansum(fij);
       sum_vs(i, :) = sum(vs(in_range, :)) ./ v0;
       Nis(i) = sum(in_range);
    end
    
    dissipations = normer(randn(N, 2)) .* [Nis, Nis] .* eta;
    repulsions = -beta.*fijtot;
    
    fs = normer(sum_vs + dissipations + repulsions) .* v0;
end


function v = normer(input)
% Returns the vectors in matrix v, normalized along the second axis.
    normval = sqrt(sum(input.^2, 2));
    v = input ./ [normval, normval];
end