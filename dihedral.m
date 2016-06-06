function phi = dihedral(xyz)
% DIHEDRAL calculate the dihedral angle between four points.
%
% phi = DIHEDRAL(XYZ) Calculate the dihedral angle for a matrix XYZ of 
% shape [3,4].
    assert(all(size(xyz) == [3,4]));
    dxyz = diff(xyz');
    dnorm = sqrt(sum(dxyz.^2, 2));
    b1 = dxyz(1,:) / dnorm(1);
    b2 = dxyz(2,:) / dnorm(2);
    b3 = dxyz(3,:) / dnorm(3);

    b12 = cross(b1, b2);
    b23 = cross(b2, b3);
    phi = atan2(cross(b12, b23) * b2', b12 * b23');
end

