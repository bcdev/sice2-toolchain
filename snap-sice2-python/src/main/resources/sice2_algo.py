import numpy as np

import sice2_constants
from sice2_constants import w, bai, xa, ya, f0, f1, f2, bet, gam, coef1, coef2, coef3, coef4
from sice2_utils import qsimp, funp


class Sice2Algo:

    def __init__(self, low_ndsi_threshold, high_ndsi_threshold):
        self.low_threshold = low_ndsi_threshold
        self.high_threshold = high_ndsi_threshold

    def compute_ndsi(self, lower_data, upper_data):
        ndsi = (upper_data - lower_data) / (upper_data + lower_data)
        return ndsi

    def compute_flags(self, ndsi):
        ndsi_low = ndsi < self.low_threshold
        ndsi_high = ndsi > self.high_threshold
        ndsi_neg = ndsi < 0.0
        ndsi_pos = ndsi >= 0.0
        ndsi_flags = ((ndsi_pos << 3) + (ndsi_neg << 2) + (ndsi_high << 1) + ndsi_low).astype(np.uint8)
        return ndsi_flags

    def ozone_scattering(self, ozone, tozon, sza, vza, toa):
        scale = np.arccos(-1.) / 180.  # rad per degree
        eps = 1.55
        # ecmwf ozone from OLCI file (in Kg.m-2) to DOBSON UNITS
        # 1 kg O3 / m2 = 46696.24  DOBSON Unit (DU)
        totadu = 46729. * ozone

        amf = 1. / np.cos(sza * scale) + 1. / np.cos(vza * scale)

        bx = (toa[20, :, :] ** (1. - eps)) * (toa[16, :, :] ** eps) / toa[6, :, :]
        o3_sice = np.log(bx) / 1.11e-4 / amf
        o3_sice[o3_sice > 500] = 999
        o3_sice[o3_sice < 0] = 999

        # Correcting TOA reflectance for ozone and water scattering

        # bav 09-02-2020: now water scattering not accounted for
        # kg/m**2. transfer to mol/cm**2
        #    roznov = 2.99236e-22  # 1 moles Ozone = 47.9982 grams
        # water vapor optical depth
        #    vap = water/roznov
        #    AKOWAT = vap/3.847e+22#    tvoda = np.exp(amf*voda*AKOWAT)

        tvoda = tozon * 0 + 1
        toa_cor_o3 = toa * np.nan

        for i in range(sice2_constants.OLCI_NUM_SPECTRAL_BANDS):
            toa_cor_o3[i, :, :] = toa[i, :, :] * tvoda[i] * np.exp(amf * tozon[i] * totadu / 404.59)

        return o3_sice, toa_cor_o3

    def view_geometry(self, vaa, saa, sza, vza):
        """
        Transfer of OLCI relative azimuthal angle to the definition used in radiative transfer code
        :param vaa:
        :param saa:
        :param sza:
        :param vza:
        :return: raa_deg, cos_sza, cos_vza, ak1, ak2, amf, raa_rt_def
        """

        raa_deg = 180. - (vaa - saa)
        sin_sza = np.sin(sza * np.pi / 180.)
        sin_vza = np.sin(vza * np.pi / 180.)

        cos_sza = np.cos(sza * np.pi / 180.)
        cos_vza = np.cos(vza * np.pi / 180.)

        ak1 = 3. * (1. + 2. * cos_sza) / 7.
        ak2 = 3. * (1. + 2. * cos_vza) / 7.

        cos_raa = np.cos(raa_deg * np.pi / 180.)
        amf = 1. / cos_sza + 1. / cos_vza
        raa_rt_def = -cos_sza * cos_vza + sin_sza * sin_vza * cos_raa

        return raa_deg, cos_sza, cos_vza, ak1, ak2, amf, raa_rt_def

    def aerosol_properties(self, aot, height, co):
        # Atmospheric optical thickness
        tauaer = aot * (w / 0.5) ** (-1.3)

        ad = height / 7400.
        ak = height * 0 + 1
        ak[ad > 1.e-6] = np.exp(-ad[ad > 1.e-6])

        taumol = np.tile(height * np.nan, (21, 1, 1))
        tau = np.tile(height * np.nan, (21, 1, 1))
        g = np.tile(height * np.nan, (21, 1, 1))
        pa = np.tile(height * np.nan, (21, 1, 1))
        p = np.tile(height * np.nan, (21, 1, 1))
        g0 = 0.5263
        g1 = 0.4627
        wave0 = 0.4685
        gaer = g0 + g1 * np.exp(-w / wave0)
        pr = 0.75 * (1. + co ** 2)

        for i in range(21):
            taumol[i, :, :] = ak * 0.00877 / w[i] ** (4.05)
            tau[i, :, :] = tauaer[i] + taumol[i, :, :]

            # aerosol asymmetry parameter
            g[i, :, :] = tauaer[i] * gaer[i] / tau[i, :, :]

            # HG phase function for aerosol
            pa[i, :, :] = (1 - g[i, :, :] ** 2) \
                          / (1. - 2. * g[i, :, :] * co + g[i, :, :] ** 2) ** 1.5

            p[i, :, :] = (taumol[i, :, :] * pr + tauaer[i] * pa[i, :, :]) / tau[i, :, :]

        return tau, p, g, gaer, taumol, tauaer

    def snow_properties(self, toa, ak1, ak2):
        # retrieval of snow properties ( R_0, size of grains from OLCI channels 865[17] and 1020nm[21]
        # assumed not influenced by atmospheric scattering and absorption processes)

        akap2 = 2.25e-6
        alpha2 = 4. * np.pi * akap2 / 1.020
        eps = 1.549559365010611

        print('rtoa[16, 350,800]: ' + str(toa[16, 350,800]))
        print('rtoa[20, 350,800]: ' + str(toa[16, 350,800]))
        print('ak1[500,500]: ' + str(ak1[500,500]))
        print('ak2[500,500]: ' + str(ak2[500,500]))
        # reflectivity of nonabsorbing snow layer
        rr1 = toa[16, :, :]
        rr2 = toa[20, :, :]
        r0 = (rr1 ** eps) * (rr2 ** (1. - eps))
        print('r0[350,800]: ' + str(r0[350,800]))

        # effective absorption length(mm)
        bal = np.log(rr2 / r0) * np.log(rr2 / r0) / alpha2 / (ak1 * ak2 / r0) ** 2
        al = bal / 1000.

        # effective grain size(mm):diameter
        D = al / 16.36
        # snow specific area ( dimension: m*m/kg)
        area = 6. / D / 0.917

        return D, area, al, r0, bal

    def prepare_coef(self, tau, g, p, am1, am2, amf, gaer, taumol, tauaer):
        astra = tau * np.nan
        rms = tau * np.nan
        t1 = tau * np.nan
        t2 = tau * np.nan

        # SOBOLEV
        oskar = 4. + 3. * (1. - g) * tau
        b1 = 1. + 1.5 * am1 + (1. - 1.5 * am1) * np.exp(-tau / am1)
        b2 = 1. + 1.5 * am2 + (1. - 1.5 * am2) * np.exp(-tau / am2)

        wa1 = 1.10363
        wa2 = -6.70122
        wx0 = 2.19777
        wdx = 0.51656
        bex = np.exp((g - wx0) / wdx)
        sssss = (wa1 - wa2) / (1. + bex) + wa2

        for i in range(21):
            astra[i, :, :] = (1. - np.exp(-tau[i, :, :] * amf)) / (am1 + am2) / 4.
            rms[i, :, :] = 1. - b1[i, :, :] * b2[i, :, :] / oskar[i, :, :] \
                           + (3. * (1. + g[i, :, :]) * am1 * am2 - 2. * (am1 + am2)) * astra[i, :, :]
            # backscattering fraction
            # t1[i, :, :] = np.exp(-(1. - g[i, :, :]) * tau[i, :, :] / am1 / 2.)
            # t2[i, :, :] = np.exp(-(1. - g[i, :, :]) * tau[i, :, :] / am2 / 2.)
            t1[i, :, :] = np.exp(-(1. - g[i, :, :]) * tau[i, :, :] / am1 / 2.
                                 / sssss[i, :, :])
            t2[i, :, :] = np.exp(-(1. - g[i, :, :]) * tau[i, :, :] / am2 / 2.
                                 / sssss[i, :, :])

        rss = p * astra
        r = rss + rms

        # SALBED
        # ratm = salbed(tau, g)
        a_s = (.18016, -0.18229, 0.15535, -0.14223)
        bs = (.58331, -0.50662, -0.09012, 0.0207)
        cs = (0.21475, -0.1, 0.13639, -0.21948)
        als = (0.16775, -0.06969, 0.08093, -0.08903)
        bets = (1.09188, 0.08994, 0.49647, -0.75218)

        a_cst = a_s[0] * g ** 0 + a_s[1] * g ** 1 + a_s[2] * g ** 2 + a_s[3] * g ** 3
        b_cst = bs[0] * g ** 0 + bs[1] * g ** 1 + bs[2] * g ** 2 + bs[3] * g ** 3
        c_cst = cs[0] * g ** 0 + cs[1] * g ** 1 + cs[2] * g ** 2 + cs[3] * g ** 3
        al_cst = als[0] * g ** 0 + als[1] * g ** 1 + als[2] * g ** 2 + als[3] * g ** 3
        bet_cst = bets[0] * g ** 0 + bets[1] * g ** 1 + bets[2] * g ** 2 + bets[3] * g ** 3

        ratm = tau * (a_cst * np.exp(-tau / al_cst) + b_cst * np.exp(-tau / bet_cst)
                      + c_cst)

        return t1, t2, ratm, r, astra, rms

    def snow_impurities(self, alb_sph, bal):
        # analysis of snow impurities
        # ( the concentrations below 0.0001 are not reliable )
        # bf    normalized absorption coefficient of pollutants ay 1000nm ( in inverse mm)
        # bm    Angstroem absorption coefficient of pollutants ( around 1 - for soot, 3-7 for dust)
        bm = np.nan * bal
        bf = bm
        p1 = bm
        p2 = bm

        ind_nonan = np.logical_and(np.logical_not(np.isnan(alb_sph[0, :, :])),
                                   np.logical_not(np.isnan(alb_sph[1, :, :])))

        p1[ind_nonan] = np.log(alb_sph[0, ind_nonan]) * np.log(alb_sph[0, ind_nonan])
        p2[ind_nonan] = np.log(alb_sph[1, ind_nonan]) * np.log(alb_sph[1, ind_nonan])
        bm[ind_nonan] = np.log(p1[ind_nonan] / p2[ind_nonan]) / np.log(w[1] / w[0])

        # type of pollutants
        ntype = np.nan * bal
        ntype[bm <= 1.2] = 1  # soot
        ntype[bm > 1.2] = 2  # dust

        soda = bm * np.nan
        soda[bm >= 0.1] = (w[0]) ** bm[bm >= 0.1]
        bf = soda * p1 / bal

        # normalized absorption coefficient of pollutants at the wavelength  1000nm
        bff = p1 / bal
        # bal   -effective absorption length in microns

        BBBB = 1.6  # enhancement factors for soot
        FFFF = 0.9  # enhancement factors for ice grains
        alfa = 4. * np.pi * 0.47 / w[0]  # bulk soot absorption coefficient at 1000nm
        DUST = 0.01  # volumetric absorption coefficient of dust

        conc = bal * np.nan
        conc[ntype == 1] = BBBB * bff[ntype == 1] / FFFF / alfa
        conc[ntype == 2] = BBBB * bff[ntype == 2] / DUST
        ntype[bm <= 0.5] = 3  # type is other or mixture
        ntype[bm >= 10.] = 4  # type is other or mixture

        return ntype, bf, conc

    def alb2rtoa(self, a, t1, t2, r0, ak1, ak2, ratm, r):
        # Function that calculates the theoretical reflectance from a snow spherical albedo a
        # This function can then be solved to find optimal snow albedo
        # Inputs:
        # a                     Surface albedo
        # r0                    reflectance of a semi-infinite non-absorbing snow layer
        #
        # Outputs:
        # rs                  surface reflectance at specific channel
        surf = t1 * t2 * r0 * a ** (ak1 * ak2 / r0) / (1 - a * ratm)
        rs = r + surf

        return rs

    def salbed(self, tau, g):
        # WARNING: NOT USED ANYMORE
        # SPHERICAL ALBEDO OF TERRESTRIAL ATMOSPHERE:
        # bav: replaced as by a_s
        # inputs:
        # tau               directional albedo ?
        # g                 asymetry coefficient
        # outputs:
        # salbed            spherical albedo
        a_s = (.18016, -0.18229, 0.15535, -0.14223)
        bs = (.58331, -0.50662, -0.09012, 0.0207)
        cs = (0.21475, -0.1, 0.13639, -0.21948)
        als = (0.16775, -0.06969, 0.08093, -0.08903)
        bets = (1.09188, 0.08994, 0.49647, -0.75218)

        a = a_s[0] * g ** 0 + a_s[1] * g ** 1 + a_s[2] * g ** 2 + a_s[3] * g ** 3
        b = bs[0] * g ** 0 + bs[1] * g ** 1 + bs[2] * g ** 2 + bs[3] * g ** 3
        c = cs[0] * g ** 0 + cs[1] * g ** 1 + cs[2] * g ** 2 + cs[3] * g ** 3
        al = als[0] * g ** 0 + als[1] * g ** 1 + als[2] * g ** 2 + als[3] * g ** 3
        bet = bets[0] * g ** 0 + bets[1] * g ** 1 + bets[2] * g ** 2 + bets[3] * g ** 3

        salbed = tau * (a * np.exp(-tau / al) + b * np.exp(-tau / bet) + c)

        return salbed

    def zbrent(self, f, x0, x1, max_iter=100, tolerance=1e-6):
        # Equation solver using Brent's method
        # https://en.wikipedia.org/wiki/Brent%27s_method
        # Brent’s is essentially the Bisection method augmented with Inverse
        # Quadratic Interpolation whenever such a step is safe. At it’s worst case
        # it converges linearly and equal to Bisection, but in general it performs
        # superlinearly; it combines the robustness of Bisection with the speedy
        # convergence and inexpensive computation of Quasi-Newtonian methods.
        # Because of this, you’re likely to find Brent’s as a default root-finding
        # algorithm in popular libraries. For example, MATLAB’s fzero, used to find
        # the root of a nonlinear function, employs a variation of Brent’s.
        # Python script from https://nickcdryan.com/2017/09/13/root-finding-algorithms-in-python-line-search-bisection-secant-newton-raphson-boydens-inverse-quadratic-interpolation-brents/

        fx0 = f(x0)
        fx1 = f(x1)

        # print(str(fx0) + ", " + str(fx1))
        if ((fx0 * fx1) > 0):
            # print("Root not bracketed "+str(fx0)+", "+str(fx1))
            # assert ((fx0 * fx1) <= 0), ("-----Root not bracketed"+str(fx0)+", "+str(fx1))
            return -999

        if abs(fx0) < abs(fx1):
            x0, x1 = x1, x0
            fx0, fx1 = fx1, fx0

        x2, fx2 = x0, fx0

        mflag = True
        steps_taken = 0
        d = np.nan

        while steps_taken < max_iter and abs(x1 - x0) > tolerance:
            fx0 = f(x0)
            fx1 = f(x1)
            fx2 = f(x2)

            if fx0 != fx2 and fx1 != fx2:
                L0 = (x0 * fx1 * fx2) / ((fx0 - fx1) * (fx0 - fx2))
                L1 = (x1 * fx0 * fx2) / ((fx1 - fx0) * (fx1 - fx2))
                L2 = (x2 * fx1 * fx0) / ((fx2 - fx0) * (fx2 - fx1))
                new = L0 + L1 + L2

            else:
                new = x1 - ((fx1 * (x1 - x0)) / (fx1 - fx0))

            if ((new < ((3 * x0 + x1) / 4) or new > x1)
                    or (mflag and (abs(new - x1)) >= (abs(x1 - x2) / 2))
                    or (mflag == False and (abs(new - x1)) >= (abs(x2 - d) / 2))
                    or (mflag and (abs(x1 - x2)) < tolerance)
                    or (mflag == False and (abs(x2 - d)) < tolerance)):
                new = (x0 + x1) / 2
                mflag = True

            else:
                mflag = False

            fnew = f(new)
            d, x2 = x2, x1

            if (fx0 * fnew) < 0:
                x1 = new
            else:
                x0 = new

            if abs(fx0) < abs(fx1):
                x0, x1 = x1, x0

            steps_taken += 1

        return x1

    def plane_albedo_sw_approx(self, D, am1):
        anka = 0.7389 - 0.1783 * am1 + 0.0484 * am1 ** 2.
        banka = 0.0853 + 0.0414 * am1 - 0.0127 * am1 ** 2.
        canka = 0.1384 + 0.0762 * am1 - 0.0268 * am1 ** 2.
        diam1 = 187.89 - 69.2636 * am1 + 40.4821 * am1 ** 2.
        diam2 = 2687.25 - 405.09 * am1 + 94.5 * am1 ** 2.
        return anka + banka * np.exp(-1000 * D / diam1) + canka \
               * np.exp(-1000 * D / diam2)

    def spher_albedo_sw_approx(self, D):
        anka = 0.6420
        banka = 0.1044
        canka = 0.1773
        diam1 = 158.62
        diam2 = 2448.18
        return anka + banka * np.exp(-1000 * D / diam1) + canka \
               * np.exp(-1000 * D / diam2)

    def qsimp(self, func, a, b):
        # integrate function between a and b using simpson's method.
        # works as fast as scipy.integrate quad
        eps = 1.e-3
        jmax = 20
        ost = -1.e30
        os = -1.e30

        for j in range(jmax):

            if (j == 0):
                st = 0.5 * (b - a) * (func(a) + func(b))
            else:
                it = 2 ** (j - 1)
                tnm = it
                delta = (b - a) / tnm
                x = a + 0.5 * delta
                sum = 0.

                for jj in range(it):
                    sum = sum + func(x)
                    x = x + delta
                st = 0.5 * (st + (b - a) * sum / tnm)
            s = (4. * st - ost) / 3.
            if (j > 4):
                if (abs(s - os) < eps * abs(os)):
                    return s
                if (s == 0) and (os == 0.):
                    return s
            os = s
            ost = st
        print("Max iteration reached")

        return s

    def BBA_calc_clean(self, al, ak1):
        # for clean snow
        # plane albedo
        sph_calc = 0  # planar

        # visible(0.3-0.7micron)

        def func_integ(x):
            return funp(x, al, sph_calc, ak1)

        p1 = qsimp(func_integ, 0.3, 0.7)

        # near-infrared (0.7-2.4micron)
        # p2 = trapzd(func_integ,0.7,2.4, 20)
        p2 = qsimp(func_integ, 0.7, 2.4)

        # spherical albedo
        sph_calc = 1  # spherical calculation

        def func_integ(x):
            return funp(x, al, sph_calc, ak1)

        # visible(0.3-0.7micron)
        # s1 = trapzd(func_integ,0.3,0.7, 20)
        s1 = qsimp(func_integ, 0.3, 0.7)
        # near-infrared (0.7-2.4micron)
        # s2 = trapzd(func_integ,0.7,2.4, 20)
        s2 = qsimp(func_integ, 0.7, 2.4)
        # shortwave(0.3-2.4 micron)
        # END of clean snow bba calculation

        return p1, p2, s1, s2

    def BBA_calc_pol(self, alb, asol, sol1_pol, sol2, sol3_pol):
        # polluted snow
        # NEW CODE FOR BBA OF BARE ICE
        # alb is either the planar or spherical albedo

        # ANAlYTICal EQUATION FOR THE NOMINATOR
        # integration over 3 segments

        # segment 1
        # QUADRATIC POLYNOMIal for the range 400-709nm
        # input wavelength
        # alam2=w[0]
        # alam3=w[5]
        # alam5=w[10]
        # alam6=w[11]
        # alam7=w[16]
        # alam8=w[20]

        alam2 = 0.4
        alam3 = 0.56
        alam5 = 0.709
        alam6 = 0.753
        alam7 = 0.865
        alam8 = 1.02

        # input reflectances
        r2 = alb[0, :]
        r3 = alb[5, :]
        r5 = alb[10, :]
        r6 = alb[11, :]
        r7 = alb[16, :]
        r8 = alb[20, :]

        # print('alb[0, 350,800]: ' + str(alb[0, 350,800]))
        # print('alb[5, 350,800]: ' + str(alb[5, 350,800]))
        # print('alb[10, 350,800]: ' + str(alb[10, 350,800]))

        sa1, a1, b1, c1 = self.quad_func(alam2, alam3, alam5, r2, r3, r5)
        ajx1 = a1 * sol1_pol
        ajx2 = b1 * coef1
        ajx3 = c1 * coef2

        aj1 = ajx1 + ajx2 + ajx3
        # segment 2.1
        # QUADRATIC POLYNOMIal for the range 709-865nm
        sa1, a2, b2, c2 = self.quad_func(alam5, alam6, alam7, r5, r6, r7)
        ajx1 = a2 * asol
        ajx2 = b2 * coef3
        ajx3 = c2 * coef4

        aj2 = ajx1 + ajx2 + ajx3  # segment 2.2
        # exponential approximation for the range 865- 2400 nm
        z1 = 0.865
        z2 = 2.4
        rati = r7 / r8
        alasta = (alam8 - alam7) / np.log(rati)
        an = 1. / alasta
        p = r7 * np.exp(alam7 / alasta)

        aj31 = (1. / an) * (np.exp(-an * z2) - np.exp(-an * z1))
        aj32 = (1. / (bet + an)) * (np.exp(-(bet + an) * z2) - np.exp(-(an + bet) * z1))
        aj33 = (1. / (gam + an)) * (np.exp(-(gam + an) * z2) - np.exp(-(an + gam) * z1))
        aj3 = (-f0 * aj31 - f1 * aj32 - f2 * aj33) * p

        BBA_vis = aj1 / sol1_pol
        BBA_nir = (aj2 + aj3) / sol2  # here segment 2.1 and 2.2 are summed
        BBA_sw = (aj1 + aj2 + aj3) / sol3_pol

        return BBA_vis, BBA_nir, BBA_sw

    def quad_func(self, x0, x1, x2, y0, y1, y2):
        # quadratic function used for the polluted snow BBA calculation
        # see BBA_calc_pol
        # compatible with arrays
        d1 = (x0 - x1) * (x0 - x2)
        d2 = (x1 - x0) * (x1 - x2)
        d3 = (x2 - x0) * (x2 - x1)

        a1 = x1 * x2 * y0 / d1 + x0 * x2 * y1 / d2 + x0 * x1 * y2 / d3
        b1 = -(x1 + x2) * y0 / d1 - (x0 + x2) * y1 / d2 - (x0 + x1) * y2 / d3
        c1 = y0 / d1 + y1 / d2 + y2 / d3
        x = x1
        sa = a1 + b1 * x + c1 * x * x

        return sa, a1, b1, c1
