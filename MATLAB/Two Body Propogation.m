%SOURCE: https://www.hlevkin.com/hlevkin/90MathPhysBioBooks/Mechanics/Curtis_OrbitamMechForEngineeringStudents.pdf

%Use Kepler's equation, solve iteratively with Newton's Method

%UNITS:
    %a: meters
    %true anomaly: radians
    %Time delta: seconds
    

function [new_semi_major_axis, new_eccentricity, new_inclination, new_ascending_node, new_periapsis, new_true_anomaly] = propogateOrbitalElements (semi_major_axis, eccentricity, inclination, ascending_node, periapsis, true_anomaly, time_delta)

%Mass of earth is extremely larger than the sattelite, so:

%gravitational  [m^3/s^2]
u = 3.986004418e14; 
a = semi_major_axis; %[m]
e = eccentricity; 
f = true_anomaly; %[radians]
dt = time_delta; %[seconds]

%convert true anomaly to eccentric anamoly
E =  2 * atan2( sqrt(1 - e) * sin(f/2), sqrt(1 + e) * cos(f/2));

%Get current true anomaly from Kepler's equation
M = E - e * sin(E);

%define mean motion
n = sqrt(u/a^3);

%Find mean anomaly at T + dt
M_new = M + n * dt;

%Use newton's method to iterate until we find the new_true_anomly. 
E_current = M_new;
tolerance = 1e-12;
MAXIMUM_ITERATIONS = 100;

for i = 1:MAXIMUM_ITERATIONS
    %Calculate the derivative of kepler's equation
    E_change = ((E_current - e * sin(E_current) - M_new)/(1 - e * cos(E_current)));
  
    %calculate next iteration
    E_current = E_current - E_change;

    %if tolerance is small enough, break out of the loop
    if(abs(E_change) < tolerance)
        break;
    end
    
end

%Nothing about the 3d Orbit changes, except the angle from the axis
new_semi_major_axis = a;
new_eccentricity = e;
new_inclination = inclination;
new_ascending_node = ascending_node;
new_periapsis = periapsis;

%updated true anomaly from our updated eccentric anomaly
new_true_anomaly = 2 * atan2( sqrt(1 + e) * sin(E_current/2), sqrt(1 - e) * cos(E_current/2));

end

%test 
[semi_major_axis, eccentricity, inclination, ascending_node, periapsis, true_anomaly] = propogateOrbitalElements((9600e3 + 21000e3)/2, 0.37255, 0, 0, 0, 0, 10800)

%correct output of 3.3371 radians


