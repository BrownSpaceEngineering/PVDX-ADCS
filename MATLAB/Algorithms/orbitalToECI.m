function [eci_x, eci_y, eci_z] = orbitalToECI(smAxis, eccentricity, inclination, aNodeLongitude, periapsisArg, trueAnomaly)
    % Converts the 6 orbital elements into ECI position vector
    % 
    % All angles must be in radians
    % Units are consistent between semi-major axis and the resultant vector
    
    % r0 is the current radius in the perifocal frame (plane of orbit)
    % since the orbit is not purely circular, it gives the current distance
    r0 = smAxis*(1-eccentricity^2)/(1+(eccentricity*cos(trueAnomaly)));
    
    % r1 is the position in the perifocal frame
    % Adjusted for the offset of the true anomaly
    % Axes: periapsis, 90 deg, normal
    r1 = [r0*cos(trueAnomaly); r0*sin(trueAnomaly); 0];

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
    
    % Rotate perifocal position with final quaternion to ECI
    [eci_x, eci_y, eci_z] = rotateVectorByQuat(r1, qComb);
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

% Test orbital
%[x,y,z] = orbitalToECI(7000, 0, 1, 0, 0, pi/2);