%INPUT:

%M_t: Magnetic field vector measured in body frame at T
%M_t_minus_1: Magnetic field vector measured in body frame at T - dT
%k: the gain value
%dT: time step

%OUPUT:
%magnetic moment to produce in body frame


function  [dipole_x, dipole_y, dipole_z] = Bdot(M_t, M_t_minus_1, k, dT)
%calculate derivative of magnetic field
bDot = (M_t - M_t_minus_1) /  dT;

%use this to find magnetic moment
m = -k * bDot;

%get vector components of magnetic moments
dipole_x = m(1);
dipole_y = m(2);
dipole_z = m(3);

end
