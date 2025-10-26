% Last edited 10/26/25 11:55 AM
function [eci_x, eci_y, eci_z, eci_vx, eci_vy, eci_vz] = orbitalToECI(smAxis, eccentricity, inclination, aNodeLongitude, periapsisArg, trueAnomaly)
    % Converts the 6 orbital elements into ECI position vector
    % 
    % All angles must be in radians
    % Units are consistent between semi-major axis and the resultant vector
    % Now position in km, velocity in km/s
    
    % r0 is the current radius in the perifocal frame (plane of orbit)
    % since the orbit is not purely circular, it gives the current distance
    r0 = smAxis*(1-eccentricity^2)/(1+(eccentricity*cos(trueAnomaly)));
    
    % r1 is the position in the perifocal frame
    % Adjusted for the offset of the true anomaly
    % Axes: periapsis, 90 deg, normal
    r1 = [r0*cos(trueAnomaly); r0*sin(trueAnomaly); 0];
    
    % Earth's gravitational parameter [km^3/s^2]
    gravConst = 398600.4418; 
    % specific angular momentum
    h = sqrt(gravConst * smAxis * (1 - eccentricity^2));
    
    % v1 is the velocity in the perifocal frame
    v1 = [-gravConst/h * sin(trueAnomaly);
           gravConst/h * (eccentricity + cos(trueAnomaly));
           0];

    % Setup for unit quaternions z and x (to rotate around axis)
    qZ = @(theta) [cos(theta/2); 0; 0; sin(theta/2)];
    qX = @(theta) [cos(theta/2); sin(theta/2); 0; 0];

    % qPeriapsisArg to orientate the x-axis at the orbital angle
    qPeriapsisArg = qZ(periapsisArg);
    % qInclination to tilt the plane relative to the equator
    qInclination = qX(inclination);
    % qLongitude to adjust to the ascending node
    qLongitude = qZ(aNodeLongitude);
    
    % Combine quaternions (Hamilton convention)
    qComb = qMult(qLongitude, qMult(qInclination, qPeriapsisArg));
    
    % Rotate perifocal position and velocity with final quaternion to ECI
    [eci_x, eci_y, eci_z] = rotateVectorByQuat(r1, qComb);
    [eci_vx, eci_vy, eci_vz] = rotateVectorByQuat(v1, qComb);
end

function qOutput = qMult(q1, q2)
    % Hamilton product of two quaternions
    w1 = q1(1); x1 = q1(2); y1 = q1(3); z1 = q1(4);
    w2 = q2(1); x2 = q2(2); y2 = q2(3); z2 = q2(4);
    
    w = w1*w2 - x1*x2 - y1*y2 - z1*z2;
    x = w1*x2 + x1*w2 + y1*z2 - z1*y2;
    y = w1*y2 - x1*z2 + y1*w2 + z1*x2;
    z = w1*z2 + x1*y2 - y1*x2 + z1*w2;
    
    qOutput = [w; x; y; z];
end

function [vX, vY, vZ] = rotateVectorByQuat(v, q)
    % Rotate vector using unit quaternion

    % Provides conjugate of q (needed for rotation)
    qConj = [q(1); -q(2:4)];
    
    % Converts the vector into a (pure) quaternion
    vQ = [0; v(:)];
    
    % v' = qvq*
    vRotQ = qMult(qMult(q, vQ), qConj);
    
    % Extracts our needed components (w should be 0 anyways)
    [vX, vY, vZ] = deal(vRotQ(2), vRotQ(3), vRotQ(4));
end

% Example Test:
% a_km = 7000;
% ecc = 0.01;
% inc_rad = 1.1;
% raan_rad = 0.3;
% argp_rad = 0.3;
% nu_rad = 0.7;
%
% Convert for keplerian2ijk 
% a_m = a_km * 1000;                  % km -> m
% inc_deg = rad2deg(inc_rad);         % rad -> deg
% raan_deg = rad2deg(raan_rad);
% argp_deg = rad2deg(argp_rad);
% nu_deg = rad2deg(nu_rad);
%
% Compute ECI using keplerian2ijk
% [rECI1_m, vECI1_ms] = keplerian2ijk(a_m, ecc, inc_deg, raan_deg, argp_deg, nu_deg);
%
% Convert keplerian2ijk outputs to km and km/s for comparison
% rECI1 = rECI1_m / 1000;      % m -> km
% vECI1 = vECI1_ms / 1000;     % m/s -> km/s
%
% Compute ECI using the quaternion script
% [x, y, z, vx, vy, vz] = orbitalToECIScript(a_km, ecc, inc_rad, raan_rad, argp_rad, nu_rad);
% rECI2 = [x, y, z];           % km
% vECI2 = [vx, vy, vz];        % km/s
%
% Print components individually 
% fprintf('rECI1: x = %.6f km, y = %.6f km, z = %.6f km\n', rECI1(1), rECI1(2), rECI1(3));
% fprintf('rECI2: x = %.6f km, y = %.6f km, z = %.6f km\n', rECI2(1), rECI2(2), rECI2(3));
%
% fprintf('vECI1: vx = %.6f km/s, vy = %.6f km/s, vz = %.6f km/s\n', vECI1(1), vECI1(2), vECI1(3));
% fprintf('vECI2: vx = %.6f km/s, vy = %.6f km/s, vz = %.6f km/s\n', vECI2(1), vECI2(2), vECI2(3));

% Print differences 
% fprintf('Position difference (ECI2 - ECI1): x = %.6f km, y = %.6f km, z = %.6f km\n', ...
%        rECI2(1)-rECI1(1), rECI2(2)-rECI1(2), rECI2(3)-rECI1(3));
% fprintf('Velocity difference (ECI2 - ECI1): vx = %.6f km/s, vy = %.6f km/s, vz = %.6f km/s\n', ...
%        vECI2(1)-vECI1(1), vECI2(2)-vECI1(2), vECI2(3)-vECI1(3));
% Example output:
% rECI1: x = 2801.905474 km, y = 3641.952685 km, z = 5209.109530 km
% rECI2: x = 2801.905474 km, y = 3641.952685 km, z = 5209.109530 km
% vECI1: vx = -6.644010 km/s, vy = -0.085065 km/s, vz = 3.698018 km/s
% vECI2: vx = -6.644010 km/s, vy = -0.085065 km/s, vz = 3.698018 km/s
% Position difference (ECI2 - ECI1): x = -0.000000 km, y = -0.000000 km, z = 0.000000 km
% Velocity difference (ECI2 - ECI1): vx = 0.000000 km/s, vy = -0.000000 km/s, vz = -0.000000 km/s