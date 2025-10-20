% Last edited 10/19/25 10:40 PM
function [sun_x, sun_y, sun_z] = sunVectorECI(currentUTC)
    % currentUTC is treated as JD
    julianDate = currentUTC;
   
    % Days since J2000.0 epoch
    julianOffset = julianDate - 2451545.0;

    % Centuries since J2000.0 epoch
    julianC = julianOffset / 36525;

    % Earth's orbital eccentricity

    % Mean anomaly and mean longitude (degrees)
    meanAnomaly = 357.529 + 35999.050*julianC;
    meanLongitude = 280.459 + 36000.770*julianC;

    % Normalize to [0, 360)
    meanAnomaly = mod(meanAnomaly, 360);
    meanLongitude = mod(meanLongitude, 360);
    
    % Sun Center Constant
    sunCenter = (1.914602 - 0.004817*julianC - 0.000014*julianC^2)*sind(meanAnomaly) + ...
    (0.019993 - 0.000101*julianC)*sind(2*meanAnomaly) + ...
    0.000289*sind(3*meanAnomaly);
    % Ecliptic Longitude (degrees)
    eclipticLongitude = meanLongitude + sunCenter;
    eclipticLongitude = mod(eclipticLongitude, 360);

    % Obliquity of ecliptic plane (degrees)
    obliquityEcliptic = 23 + 26/60 + 21.448/3600 - ...
        (46.8150*julianC + 0.00059*julianC^2 - 0.001813*julianC^3)/3600;

    % Sun direction unit vector in Earth-Centered Inertial (ECI) frame
    sun_x = cosd(eclipticLongitude);
    sun_y = cosd(obliquityEcliptic) * sind(eclipticLongitude);
    sun_z = sind(obliquityEcliptic) * sind(eclipticLongitude);
end


% Example: Current time in Julian Date

% Test of my script with current Julian Date

%jd = juliandate(datetime('now','TimeZone','UTC'));
%currentTime = jd;
%[global_x, global_y, global_z] = sunVectorECIScript(currentTime);

% Bunch of dates
%{
dates = [datetime(2020,3,15,12,0,0,'TimeZone','UTC'), ...
         datetime(2021,11,7,6,30,0,'TimeZone','UTC'), ...
         datetime(2023,5,22,18,45,0,'TimeZone','UTC'), ...
         datetime(2025,8,14,0,0,0,'TimeZone','UTC'), ...
         datetime(2029,12,31,23,59,0,'TimeZone','UTC')];
jd_array = juliandate(dates);
function dot = DotProductForDavid(x1, x2, x3, y1, y2, y3)
    dot = x1*y1 + x2*y2 + x3*y3;
end
for i = 1:length(jd_array)
    [sx, sy, sz] = sunVectorECIScript(jd_array(i));
    trueVector = planetEphemeris(jd_array(i), 'Earth', 'Sun');
    trueVector = trueVector / norm(trueVector);
    dot = DotProductForDavid(sx, sy, sz, trueVector(1), trueVector(2), trueVector(3));
    fprintf('Date: %s, Our Sun Vec: [%f, %f, %f]\n', string(dates(i)), sx, sy, sz);
    fprintf('Date: %s, MATLAB Sun Vec: [%f, %f, %f]\n', string(dates(i)), trueVector);
    fprintf('Dot: %f \n', dot);
    fprintf('________________________________________ \n');
end
%}
