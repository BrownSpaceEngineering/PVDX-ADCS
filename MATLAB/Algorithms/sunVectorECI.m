% Last edited 10/12/25 1:42 PM
function [sun_x, sun_y, sun_z] = sunVectorECI(currentUTC)
    % currentUTC is treated as JD
    julianDate = currentUTC;

    % Days since J2000.0 epoch
    julianOffset = julianDate - 2451545.0;

    % Mean anomaly and mean longitude (degrees)
    meanAnomaly = 357.529 + 0.98560028 * julianOffset;
    meanLongitude = 280.459 + 0.98564736 * julianOffset;

    % Normalize to [0, 360)
    meanAnomaly = mod(meanAnomaly, 360);
    meanLongitude = mod(meanLongitude, 360);

    % Ecliptic Longitude (degrees)
    eclipticLongitude = meanLongitude + ...
        (1.915 * sind(meanAnomaly)) + (0.020 * sind(2 * meanAnomaly));
    eclipticLongitude = mod(eclipticLongitude, 360);

    % Obliquity of ecliptic plane (degrees)
    obliquityEcliptic = 23.439 - 0.00000036 * julianOffset;

    % Sun direction unit vector in Earth-Centered Inertial (ECI) frame
    sun_x = cosd(eclipticLongitude);
    sun_y = cosd(obliquityEcliptic) * sind(eclipticLongitude);
    sun_z = sind(obliquityEcliptic) * sind(eclipticLongitude);
end


% Example: Current time in Julian Date

% currentTime = 2460634.635;

% [global_x, global_y, global_z] = sunVectorECI(currentTime);
