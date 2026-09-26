function sat = sgp4Init(oe_epoch, bstar, epochJD)
%#codegen
% sgp4Init  One-time SGP4 initialization for a new set of uplinked elements
%
% Call once each time new elements arrive (1/day), then call
% sgp4Propagate(sat, dt) at the control rate (10 Hz).
% Near-Earth branch of Vallado, Crawford, Hujsak, Kelso, "Revisiting
% Spacetrack Report #3", AIAA 2006-6753 (initl + sgp4init).
%
% INPUTS
%   oe_epoch : 6x1 SGP4/TLE *mean* elements at epoch (TEME frame)
%              [n  (rev/day)  mean motion (Kozai, as in TLE line 2)
%               e  (-)        eccentricity
%               i  (rad)      inclination
%               RAAN (rad)    right ascension of ascending node
%               argp (rad)    argument of perigee
%               M  (rad)]     mean anomaly
%   bstar    : B* drag term (1/earth radii), TLE line 1 cols 54-61
%   epochJD  : Julian date (UTC) of the element epoch. Only used for the
%              slowly varying TEME->GCRF rotation, so its precision is not
%              critical; the propagation itself only sees dt.
%
% OUTPUT
%   sat : struct of precomputed constants for sgp4Propagate.
%         sat.initErr = 0 ok, 7 = deep-space orbit (period >= 225 min),
%         not supported.

    % ================= TUNABLE CONSTANTS =================================
    % Gravity model: 72 = WGS-72 (use this for NORAD/Space-Track TLEs --
    % they are generated with WGS-72), 84 = WGS-84.
    GRAV_MODEL = 72;
    % TAI-UTC leap seconds (37 s since 2017-01-01; update if IERS adds one)
    DELTA_AT = 37;
    % =====================================================================

    if GRAV_MODEL == 84
        radiusearthkm = 6378.137;
        mu = 398600.5;
        j2 =  0.00108262998905;
        j3 = -0.00000253215306;
        j4 = -0.00000161098761;
    else % WGS-72
        radiusearthkm = 6378.135;
        mu = 398600.8;
        j2 =  0.001082616;
        j3 = -0.00000253881;
        j4 = -0.00000165597;
    end
    xke   = 60.0/sqrt(radiusearthkm^3/mu);   % sqrt(mu) in er^1.5/min
    j3oj2 = j3/j2;
    twopi = 2.0*pi;
    x2o3  = 2.0/3.0;
    temp4 = 1.5e-12;

    no_kozai = oe_epoch(1)*twopi/1440.0;     % rev/day -> rad/min
    ecco  = oe_epoch(2);
    inclo = oe_epoch(3);
    nodeo = oe_epoch(4);
    argpo = oe_epoch(5);
    mo    = oe_epoch(6);

    % ------------------------- initl ------------------------------------
    eccsq  = ecco*ecco;
    omeosq = 1.0 - eccsq;
    rteosq = sqrt(omeosq);
    cosio  = cos(inclo);
    cosio2 = cosio*cosio;

    % un-Kozai the mean motion (Brouwer mean motion)
    ak    = (xke/no_kozai)^x2o3;
    d1    = 0.75*j2*(3.0*cosio2 - 1.0)/(rteosq*omeosq);
    del   = d1/(ak*ak);
    adel  = ak*(1.0 - del*del - del*(1.0/3.0 + 134.0*del*del/81.0));
    del   = d1/(adel*adel);
    no_unkozai = no_kozai/(1.0 + del);

    ao    = (xke/no_unkozai)^x2o3;
    sinio = sin(inclo);
    po    = ao*omeosq;
    con42 = 1.0 - 5.0*cosio2;
    con41 = -con42 - cosio2 - cosio2;
    posq  = po*po;
    rp    = ao*(1.0 - ecco);

    initErr = 0;
    if twopi/no_unkozai >= 225.0
        initErr = 7;
    end

    % ------------------------- sgp4init ---------------------------------
    ss         = 78.0/radiusearthkm + 1.0;
    qzms2ttemp = (120.0 - 78.0)/radiusearthkm;
    qzms2t     = qzms2ttemp^4;

    % isimp = 1 for perigee < 220 km: drop the higher-order drag terms
    isimp = (rp < (220.0/radiusearthkm + 1.0));

    sfour  = ss;
    qzms24 = qzms2t;
    perige = (rp - 1.0)*radiusearthkm;
    if perige < 156.0
        sfour = perige - 78.0;
        if perige < 98.0
            sfour = 20.0;
        end
        qzms24 = ((120.0 - sfour)/radiusearthkm)^4;
        sfour  = sfour/radiusearthkm + 1.0;
    end
    pinvsq = 1.0/posq;

    tsi    = 1.0/(ao - sfour);
    eta    = ao*ecco*tsi;
    etasq  = eta*eta;
    eeta   = ecco*eta;
    psisq  = abs(1.0 - etasq);
    coef   = qzms24*tsi^4;
    coef1  = coef/psisq^3.5;
    cc2    = coef1*no_unkozai*(ao*(1.0 + 1.5*etasq + eeta*(4.0 + etasq)) + ...
             0.375*j2*tsi/psisq*con41*(8.0 + 3.0*etasq*(8.0 + etasq)));
    cc1    = bstar*cc2;
    cc3    = 0.0;
    if ecco > 1.0e-4
        cc3 = -2.0*coef*tsi*j3oj2*no_unkozai*sinio/ecco;
    end
    x1mth2 = 1.0 - cosio2;
    cc4    = 2.0*no_unkozai*coef1*ao*omeosq*(eta*(2.0 + 0.5*etasq) + ...
             ecco*(0.5 + 2.0*etasq) - j2*tsi/(ao*psisq)*(-3.0*con41* ...
             (1.0 - 2.0*eeta + etasq*(1.5 - 0.5*eeta)) + 0.75*x1mth2* ...
             (2.0*etasq - eeta*(1.0 + etasq))*cos(2.0*argpo)));
    cc5    = 2.0*coef1*ao*omeosq*(1.0 + 2.75*(etasq + eeta) + eeta*etasq);
    cosio4 = cosio2*cosio2;
    temp1  = 1.5*j2*pinvsq*no_unkozai;
    temp2  = 0.5*temp1*j2*pinvsq;
    temp3  = -0.46875*j4*pinvsq*pinvsq*no_unkozai;
    mdot   = no_unkozai + 0.5*temp1*rteosq*con41 + ...
             0.0625*temp2*rteosq*(13.0 - 78.0*cosio2 + 137.0*cosio4);
    argpdot = -0.5*temp1*con42 + 0.0625*temp2*(7.0 - 114.0*cosio2 + 395.0*cosio4) + ...
              temp3*(3.0 - 36.0*cosio2 + 49.0*cosio4);
    xhdot1  = -temp1*cosio;
    nodedot = xhdot1 + (0.5*temp2*(4.0 - 19.0*cosio2) + 2.0*temp3*(3.0 - 7.0*cosio2))*cosio;
    omgcof  = bstar*cc3*cos(argpo);
    xmcof   = 0.0;
    if ecco > 1.0e-4
        xmcof = -x2o3*coef*bstar/eeta;
    end
    nodecf = 3.5*omeosq*xhdot1*cc1;
    t2cof  = 1.5*cc1;
    if abs(cosio + 1.0) > 1.5e-12
        xlcof = -0.25*j3oj2*sinio*(3.0 + 5.0*cosio)/(1.0 + cosio);
    else
        xlcof = -0.25*j3oj2*sinio*(3.0 + 5.0*cosio)/temp4;
    end
    aycof  = -0.5*j3oj2*sinio;
    delmo  = (1.0 + eta*cos(mo))^3;
    sinmao = sin(mo);
    x7thm1 = 7.0*cosio2 - 1.0;

    d2 = 0.0; d3 = 0.0; d4 = 0.0;
    t3cof = 0.0; t4cof = 0.0; t5cof = 0.0;
    if ~isimp
        cc1sq = cc1*cc1;
        d2    = 4.0*ao*tsi*cc1sq;
        temp  = d2*tsi*cc1/3.0;
        d3    = (17.0*ao + sfour)*temp;
        d4    = 0.5*temp*ao*tsi*(221.0*ao + 31.0*sfour)*cc1;
        t3cof = d2 + 2.0*cc1sq;
        t4cof = 0.25*(3.0*d3 + cc1*(12.0*d2 + 10.0*cc1sq));
        t5cof = 0.2*(3.0*d4 + 12.0*cc1*d3 + 6.0*d2*d2 + 15.0*cc1sq*(2.0*d2 + cc1sq));
    end

    % TT Julian centuries since J2000 at epoch (for TEME -> GCRF)
    T0 = (epochJD + (DELTA_AT + 32.184)/86400.0 - 2451545.0)/36525.0;

    sat = struct( ...
        'initErr', initErr, 'isimp', double(isimp), ...
        'radiusearthkm', radiusearthkm, 'xke', xke, 'j2', j2, ...
        'vkmpersec', radiusearthkm*xke/60.0, ...
        'bstar', bstar, 'ecco', ecco, 'inclo', inclo, 'nodeo', nodeo, ...
        'argpo', argpo, 'mo', mo, 'no_unkozai', no_unkozai, ...
        'mdot', mdot, 'argpdot', argpdot, 'nodedot', nodedot, ...
        'nodecf', nodecf, 'cc1', cc1, 'cc4', cc4, 'cc5', cc5, ...
        't2cof', t2cof, 'omgcof', omgcof, 'xmcof', xmcof, 'eta', eta, ...
        'delmo', delmo, 'sinmao', sinmao, ...
        'd2', d2, 'd3', d3, 'd4', d4, 't3cof', t3cof, 't4cof', t4cof, 't5cof', t5cof, ...
        'con41', con41, 'x1mth2', x1mth2, 'x7thm1', x7thm1, ...
        'xlcof', xlcof, 'aycof', aycof, 'T0', T0);
end
