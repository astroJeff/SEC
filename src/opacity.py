def calc_opacity(T, rho, xmass, zion, aion, ionmax, pep_in, xne_in, eta):
    import numpy as np
    """ Calculate an approximate opacity.

    Parameters
    ----------
    T : Temperature (K)
    rho : Density (g/cm^3)
    ionmax : number of isotopes in the composition
    xmass : mass fractions of the composition
    zion : number of protons in each isotope
    aion : number of protons and neutrons in each isotope
    pep : electron-positron pressure (in erg/cm^3)
    xne : electron-positron number density (in cm^-3)
    eta : electron degeneracy parameter (chemical potential / k T)

    Returns
    -------
    orad : radiative opacity
    ocond : electron-ion opacity
    opac : total opacity in cm^2/g
    srad : radiative conductivity
    scond : conductive conductivity
    sigma : total conductivity
    """

    # constants
    con2 = 1.07726359439811217e-7
    zbound = 0.1
    t7peek = 1.0e20
    avo = 6.0221367e23
    meff = 1.194648642401440e-10
    weid = 6.884326138694269e-5
    iec = 1.754582332329132e16
    xec = 4.309054377592449e-7
    rt3 = 1.7320508075688772
    c = 2.99792458e10
    ssol = 5.67050407222e-5
    asol = 4.0*ssol/c
    k2c = (4.0/3.0)*asol*c

    tiny = 1.0e-20



    xne = xne_in + tiny
    pep = pep_in + tiny


    
    # Switch for the Iben & Christy regimes
    T6_switch1 = 1
    T6_switch2 = 1.5

    # Initialize composition variables
    w = np.zeros(6)
    zbar = 0
    ytot1 = 0
    for i in range(ionmax):
        iz = np.min([2, np.max([0, int(zion[i])])])
        ymass = xmass[i] / aion[i]
        w[iz] = w[iz] + xmass[i]
        w[iz+3] = w[iz+3] + zion[i]**2 * ymass
        zbar += zion[i] * ymass
        ytot1 += ymass

    abar = 1/ytot1
    zbar *= abar
    T6 = T/1e6
    xh = w[0]
    xhe = w[1]
    xz = w[2]




# radiative section:
# from iben apj 196 525 1975
    if (xh < 1e-5):
        xmu = max([1.0e-99, w[3]+w[4]+w[5]-1])
        xkc = (2.019e-4*rho/T6**1.7)**2.425
        xkap = 1 + xkc * (1 + xkc/24.55)
        xkb = 3.86 + 0.252*np.sqrt(xmu) + 0.018*xmu
        xka = 3.437 * (1.25 + 0.488*np.sqrt(xmu) + 0.092*xmu)
        dbar = np.exp(-xka + xkb*np.log(T6))
        oiben1 = xkap * (rho/dbar)**0.67

    if (xh<1.0e-5 or T6>T6_switch1) and (xh>=1.0e5 or xz<=zbound):
        
        if T6>T6_switch1:
            d0log = -(3.868 + 0.806*xh) + 1.8*np.log(T6)
        else:
            d0log = -(3.868 + 0.806*xh) + (3.42 - 0.52*xh)*np.log(T6)

        xka1 = 2.809 * np.exp(-(1.74 - 0.755*xh) * (np.log10(T6) - 0.22 + 0.1375*xh)**2)
        xkw  = 4.05 * np.exp(-(0.306 - 0.04125*xh) * (np.log10(T6) - 0.18 + 0.1625*xh)**2)
        xkaz = 50.0*xz*xka1 * np.exp(-0.5206*((np.log(rho)-d0log)/xkw)**2)
        dbar2log = -(4.283 + 0.7196*xh) + 3.86*np.log(T6)
        dbar1log = -5.296 + 4.833*np.log(T6)
        if (dbar2log < dbar1log): dbar1log = dbar2log
        oiben2   = (rho/np.exp(dbar1log))**(0.67) * np.exp(xkaz)


    # From Christy apj 144 108 1966
    if (T6 < T6_switch2) and (xh >= 1e-5):
        t4    = T / 1.0e4
        t4r   = np.sqrt(t4)
        t44   = t4**4
        t45   = t44 * t4
        t46   = t45 * t4
        ck1   = 2.0e6/t44 + 2.1*t46
        ck3   = 4.0e-3/t44 + 2.0e-4/rho**(0.25)
        ck2   = 4.5*t46 + 1.0/(t4*ck3)
        ck4   = 1.4e3*t4 + t46
        ck5   = 1.0e6 + 0.1*t46
        ck6   = 20.0*t4 + 5.0*t44 + t45
        xkcx  = xh*(t4r/ck1 + 1.0/ck2)
        xkcy  = xhe*(1.0/ck4 + 1.5/ck5)
        xkcz  = xz*(t4r/ck6)
        ochrs = pep * (xkcx + xkcy + xkcz)


    # Opacity in the presence of Hydrogen
    if xh >= 1e-5:
        if T6 < T6_switch1:
            orad = ochrs
        elif T6 <= T6_switch2:
            zdum  = 1.0/(T6_switch1 - T6_switch2)
            xdum  = (T6 - T6_switch2)*zdum
            ydum  = (T6 - T6_switch1)*zdum
            orad  = ochrs*xdum + oiben2*ydum
        else:
            orad = oiben2
            # Opacity in the absence of Hydrogen
    else:
       if (xz > zbound):
           orad = oiben1
       else:
           orad = oiben1*(xz/zbound) + oiben2*((zbound-xz)/zbound)


    # Add in Compton scattering opacity, Weaver et al. apj 1978 225 1021
    th = np.min([511.0, T * 8.617e-8])
    fact = 1.0 + 2.75e-2*th - 4.88e-5*th*th
    facetax = 1.0e100
    if eta<500.0:
        facetax = np.exp(0.522*eta - 1.563)
        faceta = 1.0 + facetax
        ocompt = 6.65205e-25/(fact * faceta) * xne/rho
        orad = orad + ocompt

    # Cutoff radiative opacity when 4kt/hbar is less than the plasma frequency
    T_cut = con2 * np.sqrt(xne)
    if (T < T_cut):
        if (T_cut > 200.0*T):
            orad = orad * 2.658e86
        else:
            cutfac = np.exp(T_cut/T - 1.0)
            orad = orad * cutfac

    # Fudge molecular opacity for low T
    xkf = t7peek * rho * (T * 1.0e-7)**4
    orad = xkf * orad/(xkf + orad)

    #### Conductivity section ###

    # Non-degenerate regimes, from Iben apj 196 525 1975
    log_rho = np.log10(rho)
    drel = 2.4e-7 * zbar/abar * T * np.sqrt(T)
    log_drel = np.log10(drel)
    drelim = log_drel + 1.0
    if (log_rho < drelim):
        zdel = xne/(avo*T6*np.sqrt(T6)) 
        log_zdel = np.log10(zdel)
        eta0 = np.exp(-1.20322 + (2/3) * np.log(zdel))
        eta02 = eta0*eta0

        # thpl factor
        if (log_zdel < 0.645):
            thpl = -7.5668 + np.log(zdel * (1.0 + 0.024417*zdel))
        else:
            if log_zdel < 2.5:
                thpl = -7.58110 + np.log(zdel*(1.0 + 0.02804*zdel))
                if log_zdel > 2.0:
                    thpla = thpl
                    thpl = -11.0742 + np.log(zdel**2 * (1.0 + 9.376/eta02))
                    thpl = 2.0*((2.5-log_zdel)*thpla + (log_zdel-2.0)*thpl)
            else:
                thpl = -11.0742 + np.log(zdel**2 * (1.0 + 9.376/eta02))

        # Pefac and Walf factors
        if log_zdel < 2:
            pefac = 1.0 + 0.021876*zdel
            if log_zdel > 1.5:
                pefacal = np.log(pefac)
                pefacl = np.log(0.4 * eta0 + 1.64496/eta0)
                cfac1 = 2.0 - log_zdel
                cfac2 = log_zdel - 1.5
                pefac = np.exp(2.0 * (cfac1*pefacal + cfac2*pefacl))
        else:
            pefac = 0.4 * eta0 + 1.64496/eta0

        if zdel < 40.0:
            dnefac = 1.0 + zdel * (3.4838e-4 * zdel - 2.8966e-2)
        else:
            dnefac = 1.5/eta0 * (1.0 - 0.8225/eta02)

        wpar2 = 9.24735e-3 * zdel * (rho*avo*(w[3]+w[4]+w[5])/xne + dnefac)/(np.sqrt(T6)*pefac)

        # Factor 2 error between equations a12 and a13 of Iben 1975 as
        # documented in Iben & Tutokov 370, 615, 1991, page page 621
        walf = np.log(wpar2)
        walf10 = np.log10(wpar2)


        # thx, thy, and thc factors
        if walf10 <= -3.0:
            thx = np.exp(2.413 - 0.124*walf)
        elif walf10 <= -1.0:
            thx = np.exp(0.299 - walf*(0.745 + 0.0456*walf))
        else:
            thx = np.exp(0.426 - 0.558*walf)

        if walf10 <= -3.0:
            thy = np.exp(2.158 - 0.111*walf)
        elif (walf10 <= 0.0):
            thy = np.exp(0.553 - walf*(0.55 + 0.0299*walf))
        else:
            thy = np.exp(0.553 - 0.6*walf)

        if walf10 <= -2.5:
            thc = np.exp(2.924 - 0.1*walf)
        elif walf10 <= 0.5:
            thc = np.exp(1.6740 - walf*(0.511 + 0.0338*walf))
        else:
            thc = np.exp(1.941 - 0.785*walf)

        oh = (xh*thx + xhe*thy + w[5]/3*thc) / (T6*np.exp(thpl))


    # From Yakovlev & Urpin soviet astro 1980 24 303 and
    # Potekhin et al. 1997 aa 323 415 for degenerate regimes
    if (log_rho > log_drel):
        xmas = meff * xne**(1/3)
        ymas = np.sqrt(1.0 + xmas**2)
        wfac = weid * T/ymas * xne
        cint = 1.0

        #ion-electron collision frequency and the thermal conductivity
        vie = iec * zbar * ymas * cint
        cie = wfac/vie

        # electron-electron collision frequency and thermal conductivity
        tpe = xec * np.sqrt(xne/ymas)
        yg = rt3 * tpe/T
        xrel = 1.009 * (zbar/abar * rho/1.0e6)**(1/3)
        beta2 = xrel**2/(1.0 + xrel**2)
        jy = (1 + 6/(5*xrel**2) + 2/(5*xrel**4)) * (yg**3/(3*(1+0.07414*yg)**3) * \
                                                    np.log((2.81 - 0.810*beta2 + yg)/yg) + np.pi**5/6 * yg**4/(13.91 + yg)**4 )
        vee = 0.511 * T**2 * xmas/ymas**2 * np.sqrt(xmas/ymas) * jy
        cee = wfac/vee

        # total electron thermal conductivity and conversion to an opacity
        ov1 = cie * cee/(cee + cie)
        ov = k2c/(ov1*rho) * T**3

    # Blend the opacities in the intermediate region
    if log_rho <= log_drel:
        ocond = oh
    elif (log_rho > log_drel) and (log_rho < drelim):
        farg = np.pi * (log_rho - log_drel) / 0.3
        ffac = 0.5 * (1 - np.cos(farg))
        ocond = np.exp((1-ffac)*np.log(oh) + ffac*np.log(ov))
    elif log_rho >= drelim:
        ocond = ov

    # Total opacity
    opac = orad * ocond / (ocond + orad)

    # The equivalent conductivities
    srad   = k2c/(orad*rho)  * T**3
    scond   = k2c/(ocond*rho) * T**3
    sigma   = k2c/(opac*rho)  * T**3

    return orad, ocond, opac, srad, scond, sigma
