% testSGP4  Verify sgp4Init / sgp4Propagate / sgp4Step
%
%  1) Vallado et al. 2006 verification case (sat 00005), TEME output vs
%     the published reference values (tcppver.out).
%  2) A representative 3U LEO TLE, compared against MATLAB's built-in SGP4
%     (Aerospace Toolbox propagateOrbit, PropModel="sgp4", ICRF output)
%     over 3 days -- more than the 1-day element update interval.
%  3) Round trip: osculating elements -> orbitalToECI.m -> same r, v.
%  4) Why NOT to feed outputs back in as new elements each step.
%  5) sgp4Step (Simulink wrapper) at 10 Hz, including an element update.

clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'Algorithms'));  % orbitalToECI.m

%% 1) Vallado verification case 00005 (TEME)
L1 = '1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753';
L2 = '2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667';
[oe0, bstar, epochJD] = parseTLE(L1, L2);

% tsince (min), r (km), v (km/s) from Vallado's verification output
ref = [  0.0   7022.46529266  -1400.08296755     0.03995155   1.893841015   6.405893759   4.534807250
       360.0  -7154.03120202  -3783.17682504 -3536.19412294   4.741887409  -4.151817765  -2.093935425
       720.0  -7134.59340119   6531.68641334  3260.27186483  -4.113793027  -2.911922039  -2.557327851];

fprintf('=== Test 1: Vallado verification case 00005 (TEME) ===\n');
sat = sgp4Init(oe0, bstar, epochJD);
for k = 1:size(ref,1)
    [~,~,~,~,rT,vT,err] = sgp4Propagate(sat, ref(k,1)*60);
    dr = norm(rT - ref(k,2:4)')*1e3;   % m
    dv = norm(vT - ref(k,5:7)')*1e3;   % m/s
    fprintf('  t = %6.1f min   |dr| = %.3e m   |dv| = %.3e m/s   err = %d\n', ref(k,1), dr, dv, err);
end

%% 2) 3U LEO case vs MATLAB built-in SGP4 (GCRF / ICRF)
% Representative SSO 3U cubesat TLE (~500 km, B* typical of a 3U)
L1 = tleLine('1 99999U 25001A   25266.50000000  .00005000  00000-0  30000-3 0  999');
L2 = tleLine('2 99999  97.4000 120.0000 0012000  90.0000 270.0000 15.20000000    1');
[oe0, bstar, epochJD] = parseTLE(L1, L2);

tlefile = [tempname '.tle'];
fid = fopen(tlefile, 'w'); fprintf(fid, '%s\n%s\n', L1, L2); fclose(fid);
gp = tleread(tlefile); delete(tlefile);

tMin = 0:1:3*1440;                         % 3 days at 1-min steps
jd   = epochJD + tMin/1440;
[rM, vM] = propagateOrbit(datetime(jd, 'ConvertFrom', 'juliandate'), gp, PropModel="sgp4");
rM = rM/1e3; vM = vM/1e3;                  % m -> km

sat = sgp4Init(oe0, bstar, epochJD);
rG = zeros(3, numel(jd)); vG = rG; rT = rG;
for k = 1:numel(jd)
    [rG(:,k), vG(:,k), ~, ~, rT(:,k)] = sgp4Propagate(sat, tMin(k)*60);
end
dr = vecnorm(rG - rM)*1e3;                 % m
dv = vecnorm(vG - vM)*1e3;                 % m/s
dTEME = vecnorm(rT - rM)*1e3;              % how much the frame conversion matters

fprintf('\n=== Test 2: 3U LEO vs MATLAB propagateOrbit (sgp4), 3 days ===\n');
fprintf('  GCRF position diff: max %.2f m, rms %.2f m\n', max(dr), rms(dr));
fprintf('  GCRF velocity diff: max %.4f m/s, rms %.4f m/s\n', max(dv), rms(dv));
fprintf('  (for scale: skipping TEME->GCRF would give max %.1f km error)\n', max(dTEME)/1e3);
fprintf('  Direction error of r: max %.2e deg\n', ...
    max(acosd(min(1, dot(rG, rM)./(vecnorm(rG).*vecnorm(rM))))));

figure('Name', 'SGP4 vs MATLAB');
subplot(2,1,1); plot(tMin/60, dr); grid on;
ylabel('|\Deltar| (m)'); title('sgp4Propagate vs propagateOrbit(sgp4), GCRF');
subplot(2,1,2); plot(tMin/60, dv*1e3); grid on;
ylabel('|\Deltav| (mm/s)'); xlabel('Time since epoch (h)');

%% 3) Osculating elements round trip through orbitalToECI.m
[rG1, vG1, oe] = sgp4Propagate(sat, 0.37*86400);
[x,y,z,vx,vy,vz] = orbitalToECI(oe(1), oe(2), oe(3), oe(4), oe(5), oe(6));
fprintf('\n=== Test 3: osculating elements -> orbitalToECI round trip ===\n');
fprintf('  a = %.3f km, e = %.6f, i = %.4f deg, RAAN = %.4f deg, argp = %.4f deg, nu = %.4f deg\n', ...
    oe(1), oe(2), rad2deg(oe(3)), rad2deg(oe(4)), rad2deg(oe(5)), rad2deg(oe(6)));
fprintf('  |dr| = %.3e m, |dv| = %.3e m/s\n', norm([x;y;z]-rG1)*1e3, norm([vx;vy;vz]-vG1)*1e3);

%% 4) Feeding outputs back in as the next step's elements (DON'T do this)
% Mimics a Kepler-style loop: every 0.1 s, take the elements SGP4 just
% output, treat them as a new "TLE" and propagate 0.1 s from there.
% Compared against the correct approach: fixed uplinked elements, growing dt.
dtStep = 0.1;  nSteps = 30*60/dtStep;      % 30 min at 10 Hz
mu72 = 398600.8;
oeOsc = oe0; oeMean = oe0; ep = epochJD; oscDead = 0;
for k = 1:nSteps
    if ~oscDead
        [~,~,o,~,~,~,err] = sgp4Propagate(sgp4Init(oeOsc, bstar, ep), dtStep);
        if err ~= 0 || any(~isfinite(o))
            oscDead = k;
            rOscLast = sgp4Propagate(sgp4Init(oeOsc, bstar, ep), 0);
            rOscErr  = norm(rOscLast - sgp4Propagate(sat, (k-1)*dtStep));
        else
            oeOsc = [sqrt(398600.4418/o(1)^3)*86400/(2*pi); o(2:5); o(7)];
        end
    end
    [~,~,~,m] = sgp4Propagate(sgp4Init(oeMean, bstar, ep), dtStep);
    oeMean = [sqrt(mu72/m(1)^3)*86400/(2*pi); m(2:6)];
    ep = ep + dtStep/86400;
end
rTrue = sgp4Propagate(sat, nSteps*dtStep);
rMean = sgp4Propagate(sgp4Init(oeMean, bstar, ep), 0);
fprintf('\n=== Test 4: feeding outputs back in every 0.1 s, after 30 min ===\n');
if oscDead
    fprintf('  osculating-element feedback: orbit became invalid after %.1f s (%.1f km off just before)\n', ...
        (oscDead-1)*dtStep, rOscErr);
else
    rOsc = sgp4Propagate(sgp4Init(oeOsc, bstar, ep), 0);
    fprintf('  osculating-element feedback: %.1f km error\n', norm(rOsc - rTrue));
end
fprintf('  mean-element feedback:       %.1f km error\n', norm(rMean - rTrue));

%% 5) sgp4Step at 10 Hz, with a new element set uplinked after 1 day
clear sgp4Step
t = 0:0.1:2*86400;                         % 2 days, 10 Hz
oeNew = oe0 + [1e-4; 0; 0; 0; 0; 0.5];     % stand-in for the next day's TLE
epNew = epochJD + 1;
satNew = sgp4Init(oeNew, bstar, epNew);
maxErr = 0;
for k = 1:numel(t)
    if t(k) < 86400
        [r, ~, ~, err] = sgp4Step(oe0, bstar, epochJD, t(k));
        rRef = sgp4Propagate(sat, t(k));
    else                                   % new uplink: new elements, dt restarts
        dtNew = t(k) - 86400;
        [r, ~, ~, err] = sgp4Step(oeNew, bstar, epNew, dtNew);
        rRef = sgp4Propagate(satNew, dtNew);
    end
    maxErr = max(maxErr, norm(r - rRef));
end
fprintf('\n=== Test 5: sgp4Step at 10 Hz over 2 days (%d steps) ===\n', numel(t));
fprintf('  max diff vs direct init+propagate: %.3e m, last errCode = %d\n', maxErr*1e3, err);
clear sgp4Step
t0 = tic; for k = 1:20000, sgp4Step(oe0, bstar, epochJD, k*0.1); end
fprintf('  cost per 10 Hz step (MATLAB, interpreted): %.1f us\n', toc(t0)/20000*1e6);

%% ---------------------------------------------------------------------
function [oe0, bstar, epochJD] = parseTLE(L1, L2)
% Minimal TLE parser -> inputs of sgp4Propagate
    yy  = str2double(L1(19:20));
    doy = str2double(L1(21:32));
    if yy < 57, yr = 2000 + yy; else, yr = 1900 + yy; end
    epochJD = juliandate(datetime(yr,1,1)) + doy - 1;

    mant = str2double(L1(54:59)); ex = str2double(L1(60:61));
    bstar = mant*1e-5*10^ex;

    oe0 = [str2double(L2(53:63));              % n (rev/day)
           str2double(['0.' L2(27:33)]);       % e
           deg2rad(str2double(L2(9:16)));      % i
           deg2rad(str2double(L2(18:25)));     % RAAN
           deg2rad(str2double(L2(35:42)));     % argp
           deg2rad(str2double(L2(44:51)))];    % M
end

function L = tleLine(L)
% Pad to 68 chars and append the mod-10 checksum
    L = sprintf('%-68s', L);
    c = sum(L(isstrprop(L, 'digit')) - '0') + sum(L == '-');
    L = [L num2str(mod(c, 10))];
end
