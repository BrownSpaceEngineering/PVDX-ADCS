function [r_gcrf, v_gcrf, oe_osc, oe_mean, r_teme, v_teme, errCode] = sgp4Propagate(sat, dt)
%#codegen
% sgp4Propagate  SGP4 propagation step (near-Earth / LEO only)
%
% Near-Earth branch of Vallado, Crawford, Hujsak, Kelso, "Revisiting
% Spacetrack Report #3", AIAA 2006-6753 (sgp4 routine).
%
% ALWAYS propagate from the uplinked epoch elements: dt is the total time
% since the element epoch, NOT the time since the last call. Do not feed
% the outputs back in as new elements (see testSGP4.m, test 4).
%
% INPUTS
%   sat : struct from sgp4Init (recompute only when new elements arrive)
%   dt  : time since element epoch (s), UTC-consistent clock
%
% OUTPUTS
%   r_gcrf, v_gcrf : 3x1 position (km) / velocity (km/s) in GCRF (~ICRF)
%   oe_osc  : 7x1 osculating Keplerian elements in GCRF
%             [a (km); e; i; RAAN; argp; nu (true anom); M] (rad)
%             (same convention and mu as orbitalToECI.m)
%   oe_mean : 6x1 SGP4 mean elements at dt (TEME)
%             [a (km); e; i; RAAN; argp; M] (rad) -- secular+drag only
%   r_teme, v_teme : 3x1 native SGP4 output in TEME (km, km/s)
%   errCode : 0 ok
%             1 mean e out of range   2 mean motion < 0
%             4 semi-latus rectum < 0 6 satellite has decayed
%             7 deep-space orbit (from sgp4Init), not supported

    r_gcrf = zeros(3,1); v_gcrf = zeros(3,1);
    r_teme = zeros(3,1); v_teme = zeros(3,1);
    oe_osc = zeros(7,1); oe_mean = zeros(6,1);

    if sat.initErr ~= 0
        errCode = sat.initErr;
        return;
    end

    twopi  = 2.0*pi;
    x2o3   = 2.0/3.0;
    xke    = sat.xke;
    tsince = dt/60.0;                         % minutes since epoch

    % secular gravity and atmospheric drag
    xmdf   = sat.mo + sat.mdot*tsince;
    argpdf = sat.argpo + sat.argpdot*tsince;
    nodedf = sat.nodeo + sat.nodedot*tsince;
    argpm  = argpdf;
    mm     = xmdf;
    t2     = tsince*tsince;
    nodem  = nodedf + sat.nodecf*t2;
    tempa  = 1.0 - sat.cc1*tsince;
    tempe  = sat.bstar*sat.cc4*tsince;
    templ  = sat.t2cof*t2;

    if sat.isimp == 0
        delomg   = sat.omgcof*tsince;
        delmtemp = 1.0 + sat.eta*cos(xmdf);
        delm     = sat.xmcof*(delmtemp^3 - sat.delmo);
        temp     = delomg + delm;
        mm       = xmdf + temp;
        argpm    = argpdf - temp;
        t3       = t2*tsince;
        t4       = t3*tsince;
        tempa    = tempa - sat.d2*t2 - sat.d3*t3 - sat.d4*t4;
        tempe    = tempe + sat.bstar*sat.cc5*(sin(mm) - sat.sinmao);
        templ    = templ + sat.t3cof*t3 + t4*(sat.t4cof + tsince*sat.t5cof);
    end

    nm    = sat.no_unkozai;
    em    = sat.ecco;
    inclm = sat.inclo;

    if nm <= 0.0
        errCode = 2;
        return;
    end
    am = (xke/nm)^x2o3*tempa*tempa;
    nm = xke/am^1.5;
    em = em - tempe;

    if em >= 1.0 || em < -0.001
        errCode = 1;
        return;
    end
    if em < 1.0e-6
        em = 1.0e-6;
    end
    mm    = mm + sat.no_unkozai*templ;
    xlm   = mm + argpm + nodem;

    nodem = rem(nodem, twopi);
    argpm = rem(argpm, twopi);
    xlm   = rem(xlm, twopi);
    mm    = rem(xlm - argpm - nodem, twopi);

    oe_mean = [am*sat.radiusearthkm; em; inclm; mod(nodem, twopi); ...
               mod(argpm, twopi); mod(mm, twopi)];

    sinip = sin(inclm);
    cosip = cos(inclm);

    % long period periodics (near-Earth: no lunar-solar terms)
    axnl = em*cos(argpm);
    temp = 1.0/(am*(1.0 - em*em));
    aynl = em*sin(argpm) + temp*sat.aycof;
    xl   = mm + argpm + nodem + temp*sat.xlcof*axnl;

    % solve Kepler's equation
    u    = rem(xl - nodem, twopi);
    eo1  = u;
    tem5 = 9999.9;
    ktr  = 1;
    sineo1 = 0.0; coseo1 = 1.0;
    while abs(tem5) >= 1.0e-12 && ktr <= 10
        sineo1 = sin(eo1);
        coseo1 = cos(eo1);
        tem5   = 1.0 - coseo1*axnl - sineo1*aynl;
        tem5   = (u - aynl*coseo1 + axnl*sineo1 - eo1)/tem5;
        if abs(tem5) >= 0.95
            if tem5 > 0.0
                tem5 = 0.95;
            else
                tem5 = -0.95;
            end
        end
        eo1 = eo1 + tem5;
        ktr = ktr + 1;
    end

    % short period preliminary quantities
    ecose = axnl*coseo1 + aynl*sineo1;
    esine = axnl*sineo1 - aynl*coseo1;
    el2   = axnl*axnl + aynl*aynl;
    pl    = am*(1.0 - el2);
    if pl < 0.0
        errCode = 4;
        return;
    end
    rl     = am*(1.0 - ecose);
    rdotl  = sqrt(am)*esine/rl;
    rvdotl = sqrt(pl)/rl;
    betal  = sqrt(1.0 - el2);
    temp   = esine/(1.0 + betal);
    sinu   = am/rl*(sineo1 - aynl - axnl*temp);
    cosu   = am/rl*(coseo1 - axnl + aynl*temp);
    su     = atan2(sinu, cosu);
    sin2u  = (cosu + cosu)*sinu;
    cos2u  = 1.0 - 2.0*sinu*sinu;
    temp   = 1.0/pl;
    temp1  = 0.5*sat.j2*temp;
    temp2  = temp1*temp;

    % update for short period periodics
    mrt   = rl*(1.0 - 1.5*temp2*betal*sat.con41) + 0.5*temp1*sat.x1mth2*cos2u;
    su    = su - 0.25*temp2*sat.x7thm1*sin2u;
    xnode = nodem + 1.5*temp2*cosip*sin2u;
    xinc  = inclm + 1.5*temp2*cosip*sinip*cos2u;
    mvt   = rdotl - nm*temp1*sat.x1mth2*sin2u/xke;
    rvdot = rvdotl + nm*temp1*(sat.x1mth2*cos2u + 1.5*sat.con41)/xke;

    % orientation vectors
    sinsu = sin(su);    cossu = cos(su);
    snod  = sin(xnode); cnod  = cos(xnode);
    sini  = sin(xinc);  cosi  = cos(xinc);
    xmx   = -snod*cosi;
    xmy   =  cnod*cosi;
    uvec  = [xmx*sinsu + cnod*cossu; xmy*sinsu + snod*cossu; sini*sinsu];
    vvec  = [xmx*cossu - cnod*sinsu; xmy*cossu - snod*sinsu; sini*cossu];

    r_teme = mrt*uvec*sat.radiusearthkm;
    v_teme = (mvt*uvec + rvdot*vvec)*sat.vkmpersec;

    if mrt < 1.0
        errCode = 6;
        return;
    end
    errCode = 0;

    % TEME -> GCRF, evaluated at the current time (precession ~0.14"/day)
    R = teme2gcrfMatrix(sat.T0 + dt/(86400.0*36525.0));
    r_gcrf = R*r_teme;
    v_gcrf = R*v_teme;

    oe_osc = rv2coe(r_gcrf, v_gcrf, 398600.4418);
end

% =========================================================================
function R = teme2gcrfMatrix(T)
% Rotation TEME -> GCRF using IAU-76 precession and IAU-80 nutation
% (Vallado teme2eci: r_gcrf = P * N * R3(-eqe) * r_teme). T is TT Julian
% centuries since J2000. Nutation truncated to the 10 largest IAU-80 terms
% (omitted terms each < 0.02 arcsec, < ~1 m at LEO); the ~23 mas
% FK5->GCRF frame bias is ignored.
    as2r = pi/(180.0*3600.0);
    T2 = T*T; T3 = T2*T;

    % IAU-76 precession angles
    zeta  = (2306.2181*T + 0.30188*T2 + 0.017998*T3)*as2r;
    theta = (2004.3109*T - 0.42665*T2 - 0.041833*T3)*as2r;
    z     = (2306.2181*T + 1.09468*T2 + 0.018203*T3)*as2r;

    % mean obliquity of the ecliptic
    epsb = (84381.448 - 46.8150*T - 0.00059*T2 + 0.001813*T3)*as2r;

    % IAU-80 fundamental arguments
    r  = 1296000.0;
    l  = rem(485866.733  + (1325.0*r +  715922.633)*T + 31.310*T2 + 0.064*T3, r)*as2r;
    lp = rem(1287099.804 + (  99.0*r + 1292581.224)*T -  0.577*T2 - 0.012*T3, r)*as2r;
    F  = rem(335778.877  + (1342.0*r +  295263.137)*T - 13.257*T2 + 0.011*T3, r)*as2r;
    D  = rem(1072261.307 + (1236.0*r + 1105601.328)*T -  6.891*T2 + 0.019*T3, r)*as2r;
    Om = rem(450160.280  - (   5.0*r +  482890.539)*T +  7.455*T2 + 0.008*T3, r)*as2r;

    % largest IAU-80 nutation terms, units 0.0001 arcsec
    %  [l lp F D Om], dpsi = (A + B T) sin(arg), deps = (C + D T) cos(arg)
    nut = [ 0  0  0  0  1  -171996.0 -174.2  92025.0   8.9;
            0  0  2 -2  2   -13187.0   -1.6   5736.0  -3.1;
            0  0  2  0  2    -2274.0   -0.2    977.0  -0.5;
            0  0  0  0  2     2062.0    0.2   -895.0   0.5;
            0  1  0  0  0     1426.0   -3.4     54.0  -0.1;
            1  0  0  0  0      712.0    0.1     -7.0   0.0;
            0  1  2 -2  2     -517.0    1.2    224.0  -0.6;
            0  0  2  0  1     -386.0   -0.4    200.0   0.0;
            1  0  2  0  2     -301.0    0.0    129.0  -0.1;
            0 -1  2 -2  2      217.0   -0.5    -95.0   0.3];
    dpsi = 0.0; deps = 0.0;
    for k = 1:size(nut,1)
        arg  = nut(k,1)*l + nut(k,2)*lp + nut(k,3)*F + nut(k,4)*D + nut(k,5)*Om;
        dpsi = dpsi + (nut(k,6) + nut(k,7)*T)*sin(arg);
        deps = deps + (nut(k,8) + nut(k,9)*T)*cos(arg);
    end
    dpsi = dpsi*1.0e-4*as2r;
    deps = deps*1.0e-4*as2r;
    epst = epsb + deps;

    % equation of the equinoxes (TEME -> TOD is a rotation about z)
    eqe = dpsi*cos(epsb);

    % Rk(a) rotates the coordinate frame by +a about axis k
    %  TOD   = R3(-eqe) * TEME
    %  MOD   = N * TOD,  N = R1(-epsb) * R3(dpsi) * R1(epst)
    %  J2000 = P * MOD,  P = R3(zeta) * R2(-theta) * R3(z)
    N = rot1(-epsb)*rot3(dpsi)*rot1(epst);
    P = rot3(zeta)*rot2(-theta)*rot3(z);
    R = P*N*rot3(-eqe);
end

function M = rot1(a)
    c = cos(a); s = sin(a);
    M = [1 0 0; 0 c s; 0 -s c];
end
function M = rot2(a)
    c = cos(a); s = sin(a);
    M = [c 0 -s; 0 1 0; s 0 c];
end
function M = rot3(a)
    c = cos(a); s = sin(a);
    M = [c s 0; -s c 0; 0 0 1];
end

% =========================================================================
function oe = rv2coe(r, v, mu)
% Osculating classical elements from r (km), v (km/s).
% Returns [a; e; i; RAAN; argp; nu; M] (km, rad). For near-circular or
% near-equatorial orbits the undefined angles are set to 0 and absorbed
% into nu (so orbitalToECI(a,e,i,RAAN,argp,nu) still reproduces r, v).
    small = 1.0e-10;
    rmag = norm(r); vmag = norm(v);
    h    = cross(r, v); hmag = norm(h);
    nvec = [-h(2); h(1); 0.0]; nmag = norm(nvec);
    evec = ((vmag*vmag - mu/rmag)*r - dot(r, v)*v)/mu;
    e    = norm(evec);
    xi   = 0.5*vmag*vmag - mu/rmag;
    a    = -mu/(2.0*xi);
    incl = acos(max(-1.0, min(1.0, h(3)/hmag)));

    equatorial = nmag < small*hmag;
    circular   = e < small;

    if equatorial
        RAAN = 0.0;
    else
        RAAN = atan2(nvec(2), nvec(1));
    end

    if circular
        argp = 0.0;
        if equatorial
            nu = atan2(r(2), r(1));          % true longitude
            if h(3) < 0, nu = -nu; end
        else
            nu = atan2(dot(cross(nvec, r), h)/hmag, dot(nvec, r));  % arg of latitude
        end
    else
        if equatorial
            argp = atan2(evec(2), evec(1));  % longitude of periapsis
            if h(3) < 0, argp = -argp; end
        else
            argp = atan2(dot(cross(nvec, evec), h)/hmag, dot(nvec, evec));
        end
        nu = atan2(dot(cross(evec, r), h)/hmag, dot(evec, r));
    end

    E = 2.0*atan2(sqrt(1.0 - e)*sin(nu/2.0), sqrt(1.0 + e)*cos(nu/2.0));
    M = E - e*sin(E);

    twopi = 2.0*pi;
    oe = [a; e; incl; mod(RAAN, twopi); mod(argp, twopi); mod(nu, twopi); mod(M, twopi)];
end
