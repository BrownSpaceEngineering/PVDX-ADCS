k=1; %%gain vaue?

function  [dipole_x, dipole_y, dipole_z] = Bdot(M_t, M_t_minus_1)
%calculate derivative of magnetic field
bDot = M_t - M_t_minus_1;

%use this to find magnetic moment
m = -k * bDot;

%get vector components of magnetic moments
dipole_x = m(1);
dipole_y = m(2);
dipole_z = m(3);

end
