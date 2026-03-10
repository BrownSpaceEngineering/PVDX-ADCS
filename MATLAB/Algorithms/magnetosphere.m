%NOTE: I do not truly understand this algorithm, but it seems to work and
% gives accurate results with the example cases

function B_ECI = wmmECI_embedded(r_ECI, JD)
% wmmECI_embedded  Earth's magnetic field (ECI) using embedded WMM2025 coefficients
%
% Inputs:
%   r_ECI - [3x1] ECI position (m)
%   JD    - Julian Date
%
% Output:
%   B_ECI - [3x1] magnetic field in Tesla (ECI)
%
% Reference: NOAA WMM2025 (epoch 2025.0)
% --- Convert Julian date to decimal year (keep for potential use) ---
year = 2000 + (JD - 2451545.0) / 365.25;

% --- WMM2025 embedded Gauss coefficients (nT) ---
[g,h,gdot,hdot] = loadWMM2025();

% --- Extrapolate coefficients to given year (epoch 2025.0) ---
epoch = 2025.0;
decimalYear = jd2year(JD);
g = g + (decimalYear - epoch) .* gdot;
h = h + (decimalYear - epoch) .* hdot;

% --- Compute magnetic field in ECI ---
B_ECI = computeWMMfieldFromCoeffs_manual(r_ECI, g, h, JD);

end

%% ---------------------- coefficient loader (unchanged data) ----------------------
function [g, h, gdot, hdot] = loadWMM2025()
% loadWMM2025  Embedded WMM2025 coefficients
% Output arrays: g(n+1,m+1), h(n+1,m+1), gdot(n+1,m+1), hdot(n+1,m+1)
% Units: nT and nT/yr

data = [...
 1  0  -29351.8    0.0      12.0      0.0
 1  1  -1410.8   4545.4      9.7    -21.5
 2  0  -2556.6    0.0     -11.6      0.0
 2  1   2951.1  -3133.6     -5.2    -27.7
 2  2   1649.3   -815.1     -8.0    -12.1
 3  0   1361.0    0.0       -1.3      0.0
 3  1  -2404.1    -56.6     -4.2      4.0
 3  2   1243.8    237.5      0.4     -0.3
 3  3    453.6   -549.5    -15.6     -4.1
 4  0    895.0     0.0      -1.6      0.0
 4  1    799.5    278.6     -2.4     -1.1
 4  2     55.7   -133.9     -6.0      4.1
 4  3   -281.1    212.0      5.6      1.6
 4  4     12.1   -375.6     -7.0     -4.4
 5  0   -233.2     0.0       0.6      0.0
 5  1    368.9     45.4      1.4     -0.5
 5  2    187.2    220.2      0.0      2.2
 5  3   -138.7   -122.9      0.6      0.4
 5  4   -142.0     43.0      2.2      1.7
 5  5     20.9    106.1      0.9      1.9
 6  0     64.4     0.0      -0.2      0.0
 6  1     63.8    -18.4     -0.4      0.3
 6  2     76.9     16.8      0.9     -1.6
 6  3   -115.7     48.8      1.2     -0.4
 6  4    -40.9    -59.8     -0.9      0.9
 6  5     14.9     10.9      0.3      0.7
 6  6    -60.7     72.7      0.9      0.9
 7  0     79.5     0.0      -0.0      0.0
 7  1    -77.0    -48.9     -0.1      0.6
 7  2     -8.8    -14.4     -0.1      0.5
 7  3     59.3     -1.0      0.5     -0.8
 7  4     15.8     23.4     -0.1      0.0
 7  5      2.5     -7.4     -0.8     -1.0
 7  6    -11.1    -25.1     -0.8      0.6
 7  7     14.2     -2.3      0.8     -0.2
 8  0     23.2     0.0      -0.1      0.0
 8  1     10.8     7.1       0.2     -0.2
 8  2    -17.5    -12.6      0.0      0.5
 8  3      2.0     11.4      0.5     -0.4
 8  4    -21.7     -9.7     -0.1      0.4
 8  5     16.9     12.7      0.3     -0.5
 8  6     15.0      0.7      0.2     -0.6
 8  7    -16.8     -5.2     -0.0      0.3
 8  8      0.9      3.9      0.2      0.2
 9  0      4.6      0.0     -0.0      0.0
 9  1      7.8    -24.8     -0.1     -0.3
 9  2      3.0     12.2      0.1      0.3
 9  3     -0.2      8.3      0.3     -0.3
 9  4     -2.5     -3.3     -0.3      0.3
 9  5    -13.1     -5.2      0.0      0.2
 9  6      2.4      7.2      0.3     -0.1
 9  7      8.6     -0.6     -0.1     -0.2
 9  8     -8.7      0.8      0.1      0.4
 9  9    -12.9     10.0     -0.1      0.1
 10 0     -1.3      0.0      0.1      0.0
 10 1     -6.4      3.3      0.0      0.0
 10 2      0.2      0.0      0.1     -0.0
 10 3      2.0      2.4      0.1     -0.2
 10 4     -1.0      5.3     -0.0      0.1
 10 5     -0.6     -9.1     -0.3     -0.1
 10 6     -0.9      0.4      0.0      0.1
 10 7      1.5     -4.2     -0.1      0.0
 10 8      0.9     -3.8     -0.1     -0.1
 10 9     -2.7      0.9     -0.0      0.2
 10 10    -3.9     -9.1     -0.0     -0.0
 11 0      2.9      0.0      0.0      0.0
 11 1     -1.5      0.0     -0.0     -0.0
 11 2     -2.5      2.9      0.0      0.1
 11 3      2.4     -0.6      0.0     -0.0
 11 4     -0.6      0.2      0.0      0.1
 11 5     -0.1      0.5     -0.1     -0.0
 11 6     -0.6     -0.3      0.0     -0.0
 11 7     -0.1     -1.2     -0.0      0.1
 11 8      1.1     -1.7     -0.1     -0.0
 11 9     -1.0     -2.9     -0.1      0.0
 11 10    -0.2     -1.8     -0.1      0.0
 11 11     2.6     -2.3     -0.1      0.0
 12 0     -2.0      0.0      0.0      0.0
 12 1     -0.2     -1.3      0.0     -0.0
 12 2      0.3      0.7     -0.0      0.0
 12 3      1.2      1.0     -0.0     -0.1
 12 4     -1.3     -1.4     -0.0      0.1
 12 5      0.6     -0.0     -0.0      -0.0
 12 6      0.6      0.6      0.1     -0.0
 12 7      0.5     -0.1     -0.0     -0.0
 12 8     -0.1      0.8      0.0      0.0
 12 9     -0.4      0.1      0.0     -0.0
 12 10    -0.2     -1.0     -0.1     -0.0
 12 11    -1.3      0.1     -0.0      0.0
 12 12    -0.7      0.2     -0.1     -0.1 ];

maxN = max(data(:,1));
g = zeros(maxN+1);
h = zeros(maxN+1);
gdot = zeros(maxN+1);
hdot = zeros(maxN+1);

for i = 1:size(data,1)
    n = data(i,1);
    m = data(i,2);
    g(n+1,m+1) = data(i,3);
    h(n+1,m+1) = data(i,4);
    gdot(n+1,m+1) = data(i,5);
    hdot(n+1,m+1) = data(i,6);
end
end

%% ---------------------- core compute routine ----------------------
function B_ECI = computeWMMfieldFromCoeffs_manual(r_ECI, g, h, jd)
    % Convert ECI -> ECEF
    theta = gmst_from_jd(jd);
    Rz = [ cos(theta)  sin(theta) 0;
          -sin(theta)  cos(theta) 0;
             0           0        1];

    r_ecef = Rz * r_ECI(:);

    utc = datetime('now','TimeZone','UTC');
    [r_eci_test] = ecef2eci(utc,r_ecef);

    % ECEF -> geodetic (WGS84)
    [lat, lon, alt] = ecef2geodetic(r_ecef);
    ecef = lla2ecef([rad2deg(lat), rad2deg(lon), alt]);
    
    % Synthesize magnetic field in NED (nT)
    [B_north, B_east, B_down] = synthesizeMagField_manual(lat, lon, alt, g, h);

    % Test with Matlab's script
    lat_deg = rad2deg(lat); 
    lon_deg = rad2deg(lon);
    [XYZ] = wrldmagm(alt, lat_deg, lon_deg, decyear(2025,11,16),'2025');

    fprintf('Matlab NED Vector (Correct magnetic field): [%.4e, %.4e, %.4e]\n', XYZ(1), XYZ(2), XYZ(3));
    fprintf('Our NED Vector: [%.4e, %.4e, %.4e]\n', B_north, B_east, B_down);
    accuracy = dot([XYZ(1), XYZ(2), XYZ(3)]/norm([XYZ(1), XYZ(2), XYZ(3)]), [B_north, B_east, B_down]/norm([B_north, B_east, B_down]));
    angle_error_deg = rad2deg(acos(accuracy));
    fprintf('Our accuracy: %f\n', accuracy);
    fprintf('Our angular error (degrees): %f\n', angle_error_deg);
    % NED -> ECEF (NED vector to ECEF vector)
    R_ned2ecef = [ -sin(lat)*cos(lon), -sin(lat)*sin(lon),  cos(lat);
                   -sin(lon),            cos(lon),           0;
                   -cos(lat)*cos(lon),  -cos(lat)*sin(lon), -sin(lat)];
    B_ecef = R_ned2ecef * [B_north; B_east; B_down];
    
    B_ECI = Rz' * B_ecef;    % still in nT
    B_ECI = B_ECI * 1e-9;    % convert nT -> T
end

%% ---------------------- Julian -> decimal year ----------------------
function year = jd2year(jd)
    % Convert Julian Date to decimal year
    days_per_year = 365.25;
    year = 2000 + (jd - 2451545.0)/days_per_year; % JD2000 reference
end

%% ---------------------- GMST ----------------------
function gmst_rad = gmst_from_jd(JD)
    % Returns Greenwich Mean Sidereal Time (gmst) in radians
    % Implementation based on IAU conventions (approximate but accurate for our needs)
    T = (JD - 2451545.0)/36525;
    GMST_sec = 67310.54841 + (876600*3600 + 8640184.812866)*T + 0.093104*T.^2 - 6.2e-6*T.^3;
    GMST_deg = mod(GMST_sec/240, 360); 
    gmst_rad = deg2rad(GMST_deg);
end

%% ---------------------- ECEF -> geodetic (WGS-84) ----------------------
function [lat, lon, alt] = ecef2geodetic(r_ecef)
    a = 6378137.0;              % semi-major axis (m)
    f = 1/298.257223563;        % flattening
    b = a * (1 - f);            % semi-minor axis
    e2 = f*(2-f);               % first eccentricity squared
    ep2 = e2/(1-e2);            % second eccentricity squared (e'^2)
    
    x = r_ecef(1);
    y = r_ecef(2);
    z = r_ecef(3);
    
    lon = atan2(y, x);
    
    % Distance from z-axis
    p = sqrt(x^2 + y^2);
    
    % Handle pole case
    if p < 1e-12
        lon = 0;
        if z >= 0
            lat = pi/2;
            alt = z - b;
        else
            lat = -pi/2;
            alt = -z - b;
        end
        return;
    end
    
    % Bowring's initial estimate for parametric latitude
    theta = atan2(z * a, p * b);
    
    % Iterate to find geodetic latitude
    for iter = 1:3  % Usually converges by 3 iterations
        sin_theta = sin(theta);
        cos_theta = cos(theta);
        
        % Calculate geodetic latitude
        num = z + ep2 * b * sin_theta^3;
        den = p - e2 * a * cos_theta^3;
        lat = atan2(num, den);
        
        % Update theta (parametric latitude) for next iteration
        sin_lat = sin(lat);
        cos_lat = cos(lat);
        N = a / sqrt(1 - e2 * sin_lat^2);
        
        % Compute altitude (Bowring 1985 formula)
        alt = p * cos_lat + z * sin_lat - a * sqrt(1 - e2 * sin_lat^2);
        
        % New theta estimate
        theta_new = atan2(z * (1 - e2 * N / (N + alt)), p);
        
        % Check convergence
        if abs(theta_new - theta) < 1e-14
            break;
        end
        theta = theta_new;
    end
end

%% ---------------------- Synthesize magnetic field (Schmidt semi-normalized) ----------------------
function [B_N, B_E, B_D] = synthesizeMagField_manual(lat, lon, alt, g, h)
    % lat, lon in radians
    % alt in meters
    % g, h are Schmidt semi-normalized coefficients (n x m)
    
    % 1. CONSTANTS (WGS84)
    a = 6371200;        % Geomagnetic Reference Radius (m)
    b = 6356752.3142;   % Polar Radius (m)
    Re = 6378137;       % Equatorial Radius (m)
    f = 1/298.257223563; % Flattening
    
    % 2. CONVERT GEODETIC TO GEOCENTRIC
    % Calculate Geocentric Latitude (phi_prime) and Radius (r)
    sin_lat = sin(lat);
    cos_lat = cos(lat);
    rc = Re / sqrt(1 - (2*f - f*f)*sin_lat^2); % Radius of curvature
    
    p = (rc + alt) * cos_lat;
    z = (rc * (1 - f)^2 + alt) * sin_lat;
    
    r = sqrt(p^2 + z^2);           % Geocentric radial distance
    phi_prime = asin(z / r);       % Geocentric latitude
    
    % The angle difference for vector rotation later
    psi = lat - phi_prime; 
    
    % 3. SETUP SPHERICAL HARMONICS
    theta = pi/2 - phi_prime; % Colatitude (geocentric)
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    
    nmax = size(g, 1) - 1;
    P = zeros(nmax+1, nmax+1);
    dP = zeros(nmax+1, nmax+1);
    
    % 4. LEGENDRE POLYNOMIALS (Schmidt Semi-Normalized)
    P(1,1) = 1; 
    P(2,1) = cos_theta;               % P(1,0)
    P(2,2) = sin_theta;               % P(1,1)
    
    dP(1,1) = 0;
    dP(2,1) = -sin_theta;             % dP(1,0)
    dP(2,2) = cos_theta;              % dP(1,1)

    for n = 2:nmax
        % Recursive coefficients
        for m = 0:n
            ind_n = n+1; 
            ind_m = m+1;
            
            if m == n
                % Diagonal: P(n,n)
                P(ind_n, ind_m) = sin_theta * P(n, n) * sqrt(1 - 1/(2*n));
                dP(ind_n, ind_m) = sin_theta * dP(n, n) * sqrt(1 - 1/(2*n)) + ...
                                   cos_theta * P(n, n) * sqrt(1 - 1/(2*n));
            else
                % Standard Recurrence for Schmidt
                % We need P(n-1, m) and P(n-2, m)
                
                % Coefficient K (different for m=0 vs m>0 due to Schmidt)
                if m == 0
                end
                
                
                denom = sqrt(n^2 - m^2);
                K = (2*n - 1) / denom;
                M = sqrt(((n-1)^2 - m^2)) / denom;
                
                P(ind_n, ind_m) = K * cos_theta * P(n, ind_m) - M * P(n-1, ind_m);
                
                % Derivative
                dP(ind_n, ind_m) = K * (cos_theta * dP(n, ind_m) - sin_theta * P(n, ind_m)) ...
                                   - M * dP(n-1, ind_m);
            end
        end
    end
    
    % 5. SUMMATION
    Br = 0;     % -d/dr (Radial component)
    Bt = 0;     % -1/r d/dtheta (Theta component, South-ish)
    Bphi = 0;   % -1/(r sin theta) d/dphi (East component)
    
    % Precompute radii powers
    
    for n = 1:nmax
        ratio = (a / r)^(n + 2);
        
        for m = 0:n
            ind_n = n+1;
            ind_m = m+1;
            
            gm = g(ind_n, ind_m);
            hm = h(ind_n, ind_m);
            
            sin_mlon = sin(m * lon);
            cos_mlon = cos(m * lon);
            
            % Calculate gauss sums
            xy_comp = (gm * cos_mlon + hm * sin_mlon);
            z_comp  = (gm * sin_mlon - hm * cos_mlon);
            
            % Radial (Up/Down axis in Geocentric)
            Br = Br + (n + 1) * ratio * xy_comp * P(ind_n, ind_m);
            
            % Theta (South axis in Geocentric)
            Bt = Bt + ratio * xy_comp * dP(ind_n, ind_m);
            
            % Phi (East axis)
            if m > 0 && abs(sin_theta) > 1e-10
                 Bphi = Bphi + ratio * m * z_comp * P(ind_n, ind_m) * (1/sin_theta);
            end
        end
    end
    
    % Result is currently:
    % Br   (Geocentric Radial, points OUT)
    % Bt   (Geocentric Theta, points SOUTH)
    % Bphi (Geocentric Phi, points EAST)
    
    % 6. ROTATE TO GEODETIC (NED)
    % B_radial points OUT. B_Down points IN. So B_Down_Geocentric = -Br.
    % B_North_Geocentric = -Bt.
    
    B_X_geo = -Bt; % North (Geocentric)
    B_Y_geo = Bphi;% East (Geocentric)
    B_Z_geo = -Br; % Down (Geocentric)
    
    % Rotate by psi (difference between geodetic and geocentric latitude)
    % [ B_North ]   [ cos(psi)  -sin(psi) ] [ B_X_geo ]
    % [ B_Down  ]   [ sin(psi)   cos(psi) ] [ B_Z_geo ]
    
    B_N = -B_X_geo * cos(psi) - B_Z_geo * sin(psi);
    B_E = B_Y_geo; % East is invariant to latitude rotation
    B_D = B_X_geo * sin(psi) + B_Z_geo * cos(psi);
    % Return nT components
end

% --- TEST CASES ---
utcNow = datetime('now','TimeZone','UTC');
jd = juliandate(utcNow);
fprintf('Our Julian Day: [%f]\n', jd);

% Define Test Cases (All units in Meters)
tests = [
    % 1. SAA CENTER (The Danger Zone)
    % Location: ~30°S, 40°W (South Atlantic Anomaly). Weakest field, high radiation.
    struct('name', 'SAA Center', 'r', [4.5e6; -3.7e6; -3.4e6]),

    % 2. ISS ORBIT MAX LATITUDE
    % Location: 51.6°N, 0°E. The highest latitude a standard ISS-deployed CubeSat reaches.
    struct('name', 'ISS Max Lat', 'r', [4.2e6; 0; 5.3e6]),

    % 3. MAGNETIC NORTH POLE (Approx)
    % Location: High Arctic (~85°N). Field lines point straight down (-Z).
    struct('name', 'Mag N Pole', 'r', [0.5e6; 0; 6.75e6]),

    % 4. MAGNETIC SOUTH POLE (Approx)
    % Location: Off the coast of Antarctica (~65°S, 135°E). Field lines point straight up (+Z).
    struct('name', 'Mag S Pole', 'r', [-2.0e6; 2.0e6; -6.1e6]),

    % 5. THE "TERMINATOR" (Equatorial Crossing)
    % Location: 0°N, 90°E (Indian Ocean). Simple check for horizontal field lines.
    struct('name', 'Eq 90 deg E', 'r', [0; 6.77e6; 0]),

    % 6. BERMUDA (North Atlantic)
    % Location: ~30°N, 65°W. A standard mid-latitude verification point.
    struct('name', 'Bermuda',   'r', [2.5e6; -5.3e6; 3.4e6]),

    % 7. PACIFIC "VOID"
    % Location: 0°N, 180°E. The middle of the Pacific. Good for checking sign flips on Longitude.
    struct('name', 'Pacific 180', 'r', [-6.77e6; 0; 0]),

    % 8. 45-DEGREE TEST
    % Location: 45°N, 45°E. X, Y, and Z coords are roughly equal magnitude. 
    % Good for checking matrix mixing errors.
    struct('name', '45N 45E',     'r', [3.4e6; 3.4e6; 4.8e6]),

    % 9. POLAR ORBIT (Sun-Synchronous)
    % Location: 98° Inclination crossing the equator. (Actually, let's do high lat 80°S).
    % Many CubeSats use SSO orbits.
    struct('name', 'SSO Polar',   'r', [1.1e6; 0; -6.6e6]),

    % 10. VANDENBERG LAUNCH SITE
    % Location: 34°N, 120°W (California). Common launch site for polar CubeSats.
    struct('name', 'Vandenberg',  'r', [-2.8e6; -4.8e6; 3.8e6])
];

% --- EXECUTION LOOP ---
fprintf('----------------------------------------------------------------\n');
fprintf('%-15s | %-12s | %-30s\n', 'Test Name', 'Input Mag', 'Output B_ECI [Bx, By, Bz] (nT)');
fprintf('----------------------------------------------------------------\n');

for i = 1:length(tests)
    r_check = tests(i).r;
    
    % Run your function
    B_vec = wmmECI_embedded(r_check, jd);
    
    % Calculate Magnitude for quick check
    B_mag = norm(B_vec)*10^9;
    
    fprintf('Our Starting ECI position: [%.4e, %.4e, %.4e]\n', r_check);
    % Print results
    fprintf('%-15s | %6.0f km   | [%8.1f, %8.1f, %8.1f] (Total: %.1f nT)\n', ...
        tests(i).name, ...
        norm(r_check)/1000, ...
        B_vec(1)*10^9, B_vec(2)*10^9, B_vec(3)*10^9, B_mag);
    fprintf('----------------------------------------------------------------\n');

end