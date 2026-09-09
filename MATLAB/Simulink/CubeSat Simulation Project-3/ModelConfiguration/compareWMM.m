function results = compareWMM(userFcn, opts)
%compareWMM Validate a user WMM2025 implementation against MATLAB's wrldmagm.
%
%   results = compareWMM(@myWMM) compares myWMM against wrldmagm over a
%   global lat/lon/alt/time grid and reports per-component error statistics.
%
%   userFcn must have the signature:
%       B_ned_nT = userFcn(height_m, lat_deg, lon_deg, decimalYear)
%   returning a 3-element NED vector in nanotesla.  Use opts to adapt if
%   your function differs (see below) -- unit/frame mismatch is the single
%   most common cause of a false failure.
%
%   opts (name-value):
%       Scale       multiply your output by this to reach nT   (default 1)
%                   tesla -> 1e9, gauss -> 1e5, microtesla -> 1e3
%       Frame       'NED' (default) | 'ENU' | 'NEU'
%       AltUnits    'm' (default) | 'km'
%       TimeFormat  'decyear' (default) | 'juliandate' | 'datetime'
%       Grid        'coarse' (default, ~500 pts) | 'fine' (~5000 pts)
%       Tol         pass threshold on max component error, nT (default 1)
%
%   Reference: WMM2025, valid 2025.0 - 2030.0.

arguments
    userFcn (1,1) function_handle
    opts.Scale (1,1) double = 1
    opts.Frame (1,1) string {mustBeMember(opts.Frame,["NED","ENU","NEU"])} = "NED"
    opts.AltUnits (1,1) string {mustBeMember(opts.AltUnits,["m","km"])} = "m"
    opts.TimeFormat (1,1) string {mustBeMember(opts.TimeFormat,["decyear","juliandate","datetime"])} = "decyear"
    opts.Grid (1,1) string {mustBeMember(opts.Grid,["coarse","fine"])} = "coarse"
    opts.Tol (1,1) double = 1
end

% -- Build the test grid ------------------------------------------------
if opts.Grid == "coarse"
    lat = -80:20:80;  lon = -180:30:150;  alt = [0 400e3 800e3];  yr = [2025.0 2027.5 2029.9];
else
    lat = -89:5:89;   lon = -180:10:170;  alt = [0 200e3 400e3 800e3];  yr = [2025.0 2026.5 2027.5 2029.9];
end
[LA,LO,AL,YR] = ndgrid(lat,lon,alt,yr);
n = numel(LA);

ref  = zeros(n,3);
mine = zeros(n,3);
fail = strings(0,1);

for k = 1:n
    % --- reference: WMM2025 via wrldmagm (geodetic lat, ellipsoidal ht, decyear)
    ref(k,:) = wrldmagm(AL(k), LA(k), LO(k), YR(k), '2025')';

    % --- candidate, with interface adaptation
    h = AL(k);  if opts.AltUnits == "km", h = h/1e3; end
    switch opts.TimeFormat
        case "decyear",     t = YR(k);
        case "juliandate",  t = juliandate(decyear2dt(YR(k)));
        case "datetime",    t = decyear2dt(YR(k));
    end
    try
        v = userFcn(h, LA(k), LO(k), t);
    catch ME
        fail(end+1) = sprintf('lat %g lon %g alt %g yr %g: %s', LA(k),LO(k),AL(k),YR(k),ME.message); %#ok<AGROW>
        mine(k,:) = NaN;  continue
    end
    v = double(v(:))' * opts.Scale;
    switch opts.Frame                       % re-express as NED
        case "ENU", v = [v(2), v(1), -v(3)];
        case "NEU", v = [v(1), v(2), -v(3)];
    end
    mine(k,:) = v;
end

% -- Statistics ---------------------------------------------------------
err  = mine - ref;
aerr = abs(err);
results.grid   = table(LA(:),LO(:),AL(:),YR(:),'VariableNames',{'lat','lon','alt_m','decyear'});
results.ref_nT  = ref;
results.mine_nT = mine;
results.err_nT  = err;
results.errors  = fail;

comp = ["North","East","Down"];
fprintf('\n=== WMM2025 comparison: %d points ===\n', n);
fprintf('%-8s %12s %12s %12s\n','comp','max|err| nT','rms nT','max rel');
for i = 1:3
    rel = aerr(:,i) ./ max(abs(ref(:,i)),1);
    fprintf('%-8s %12.4g %12.4g %12.3e\n', comp(i), max(aerr(:,i)), rms(err(~isnan(err(:,i)),i)), max(rel));
end
magRef  = vecnorm(ref,2,2);
magMine = vecnorm(mine,2,2);
fprintf('%-8s %12.4g %12.4g %12.3e\n','|B|', max(abs(magMine-magRef)), rms(magMine-magRef,'omitnan'), max(abs(magMine-magRef)./magRef));

results.maxErr = max(aerr(:));
results.pass   = results.maxErr <= opts.Tol && isempty(fail);
if ~isempty(fail)
    fprintf('\n%d point(s) errored, first: %s\n', numel(fail), fail(1));
end
fprintf('\nVERDICT: %s  (max component error %.4g nT, tol %g)\n\n', ...
    string(results.pass).replace("true","PASS").replace("false","FAIL"), results.maxErr, opts.Tol);

% -- Worst offenders ----------------------------------------------------
if ~results.pass
    [~,ix] = maxk(max(aerr,[],2), min(5,n));
    fprintf('Worst points:\n');
    disp([results.grid(ix,:), table(max(aerr(ix,:),[],2),'VariableNames',{'maxAbsErr_nT'})]);
end
end

function dt = decyear2dt(y)
    y0 = floor(y);
    d0 = datetime(y0,1,1);
    dt = d0 + years(0) + (datetime(y0+1,1,1)-d0)*(y-y0);
end
