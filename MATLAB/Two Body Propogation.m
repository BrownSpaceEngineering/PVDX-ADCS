%SOURCE: https://www.hlevkin.com/hlevkin/90MathPhysBioBooks/Mechanics/Curtis_OrbitamMechForEngineeringStudents.pdf

%Use Kepler's equation, solve iteratively with Newton's Method

%UNITS:
    %a: meters
    %inclination, ascending_node, periapsis, true anomaly: degrees
    %Time delta: seconds

function [new_semi_major_axis, new_eccentricity, new_inclination, new_ascending_node, new_periapsis, new_true_anomaly] = propogateOrbitalElements (semi_major_axis, eccentricity, inclination, ascending_node, periapsis, true_anomaly, time_delta)

%Mass of earth is extremely larger than the sattelite, so:

%gravitational  
u = 3.986004418e14; %[m^3/s^2]
a = semi_major_axis; %[m]
e = eccentricity; 
f = deg2rad(true_anomaly); %[radians]
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
MAXIMUM_ITERATIONS = 10;

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
new_true_anomaly = rad2deg(new_true_anomaly);
end

%TESTING ON TEXTBOOK CALCULATION: https://www.hlevkin.com/hlevkin/90MathPhysBioBooks/Mechanics/Curtis_OrbitamMechForEngineeringStudents.pdf
[semi_major_axis, eccentricity, inclination, ascending_node, periapsis, true_anomaly] = propogateOrbitalElements((9600e3 + 21000e3)/2, 0.37255, 0, 0, 0, 0, 10800);
%correct output of 3.3371 radians

%TESTING ON MATLAB LIBRARY: https://www.mathworks.com/help/aerotbx/ug/propagateorbit.html

function testOnMatlab(axis, e0, inc0, RAAN0, argp0, nu0)
% 1 min increments of testing, over 2 hours
sampleTime = 60;
startTime = datetime(2025,10,22,12,0,0); 
endTime = datetime(2025,10,22,14,0,0); 
time = startTime:seconds(sampleTime):endTime;

%Testing on Matlab aerospace toolbox functions
[position, velocity] = propagateOrbit(time, axis, e0, inc0, RAAN0, argp0, nu0, PropModel="two-body-keplerian"); 
[a, ecc, incl, RAAN, argp, nu, truelon, arglat, lonper] = ijk2keplerian(position, velocity); 

true_anamoly_vector = [];
%loop over every minute, to get the next prediction
for j=0:120
    [a1, ecc1, incl1, RAAN1, argp1, nu1] = propogateOrbitalElements(axis, e0, inc0, RAAN0, argp0, nu0, j*sampleTime);
    true_anamoly_vector(j+1) = nu1;
    
end
%get total error between the verified function and mine
error = abs(wrapTo360(true_anamoly_vector) - wrapTo360(nu));
mean(error);
end

%testing various e
testOnMatlab(1e7, 0.7, 70, 10, 120, 1); %mean error e-8
testOnMatlab(1e7, 0.5, 70, 10, 120, 1); %mean error e-8
testOnMatlab(1e7, 0.10, 70, 10, 120, 1); %mean error e-9
testOnMatlab(1e7, 0.9, 70, 10, 120, 1); %mean error e-9
testOnMatlab(1e7, 0.95, 70, 10, 120, 1); %mean error e-9
testOnMatlab(1e7, 0.05, 70, 10, 120, 1); %mean error e-8

%testing various true anamolys
testOnMatlab(1e7, 0.3, 150, 10, -40, 1); %mean error e-8
testOnMatlab(1e7, 0.3, 150, 10, -40, -180); %mean error e-8
testOnMatlab(1e7, 0.3, 150, 10, -40, 90); %mean error e-8
testOnMatlab(1e7, 0.3, 150, 10, -40, 240); %mean error e-8

%testing various a
testOnMatlab(1e10, 0.7, 78, 13, 70, 28); %mean error e-9
testOnMatlab(1e3, 0.5, 70, 10, 200, 59); %mean error e-8
testOnMatlab(1e12, 0.10, 70, 10, 40, 192); %mean error e-8
testOnMatlab(1e6, 0.9, 70, 10, 10, -25); %mean error e-8
testOnMatlab(1e8, 0.95, 70, 10, 120, 320); %mean error e-8
testOnMatlab(1e7, 0.05, 70, 10, 80, 143); %mean error e-8

