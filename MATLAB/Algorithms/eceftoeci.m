% This function can actually act as our PVD down vector, if we note the
% following

% Providence's lat-lon is 41.825226, -71.418884
% Converted to ECEF, this is lla2ecef(deg2rad([41.825226, -71.418884, 0]))
% which is 6.3761e6, -0.1387e6, 0.0807e6. This is a constant that we don't
% need to recalculate


function [r_eci] = eceftoeci(r_ecef, unix_seconds)
    % Step 1: Unix time -> Julian Date
    JD = unix_to_jd(unix_seconds);

    % Step 2: Julian Date -> GMST (degrees)
    theta_deg = jd_to_gmst(JD);

    % Step 3: Build rotation matrix R_z(theta)
    theta = deg2rad(theta_deg);
    R = [ cos(theta), -sin(theta), 0; ...
          sin(theta),  cos(theta), 0; ...
          0,           0,          1];

    % Step 4: Rotate position
    r_eci = R * r_ecef;

end


% -------------------------------------------------------------------------
function JD = unix_to_jd(unix_seconds)
% UNIX_TO_JD Convert Unix timestamp to Julian Date.
%   Unix epoch (1970-01-01 00:00:00 UTC) = JD 2440587.5
    JD = unix_seconds / 86400.0 + 2440587.5;
end


% -------------------------------------------------------------------------
function theta_deg = jd_to_gmst(JD)
% JD_TO_GMST Compute Greenwich Mean Sidereal Time (degrees) from Julian Date.

    % Julian centuries since J2000.0
    T = (JD - 2451545.0) / 36525.0;

    % GMST in degrees
    theta_deg = 280.46061837 ...
              + 360.98564736629 * (JD - 2451545.0) ...
              + 0.000387933     * T^2 ...
              - T^3             / 38710000.0;

    % Reduce to [0, 360)
    theta_deg = mod(theta_deg, 360.0);
end

% Put this in a separate function to test. We get ~7km error at worst,
% which isn't awful

% function test()
%     worst_err = 0;
%     for i = 1:100
%         coords = randi([600000, 650000], 3, 1);
%         time = 1772814426;
%         our_r = eceftoeci(coords, time);
%         utc_time = datetime(time, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
%         matlab_r = ecef2eci(utc_time, coords);
%         worst_err = max(worst_err, norm((our_r - matlab_r)));
%     end
%     disp(worst_err);
% end