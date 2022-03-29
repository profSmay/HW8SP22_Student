import numpy as np
from scipy.interpolate import griddata
from scipy.optimize import fsolve
from pyXSteam.XSteam import XSteam

class UnitConverter():
    def __init__(self):
        """
        This unit converter class is useful for the pipe network and perhaps other problems.
        The strategy is (number in current units)*(conversion factor)=(number desired units), for instance:
            1(ft)*(self.ft_to_m) = 1/3.28084 (m)
            1(in^2)*(self.in2_to_m2) = 1*(1/(12*3.28084))**2 (m^2)
        """

    """
    These constants can be used directly from the class without instantiating an object
    """
    ft_to_m = 1 / 3.28084
    ft2_to_m2 = ft_to_m ** 2
    ft3_to_m3 = ft_to_m ** 3
    m3_to_ft3=1/ft3_to_m3
    ft3_to_L = ft3_to_m3 * 1000
    L_to_ft3 = 1 / ft3_to_L
    in_to_m = ft_to_m / 12
    m_to_in = 1 / in_to_m
    in2_to_m2 = in_to_m ** 2
    m2_to_in2 = 1 / in2_to_m2
    g_SI = 9.80665  # m/s^2
    g_EN = 32.174  # 32.174 ft/s^2
    gc_EN = 32.174  # lbm*ft/lbf*s^2
    gc_SI = 1.0  # kg*m/N*s^2
    lbf_to_kg = 1 / 2.20462
    kg_to_lbf = 1/lbf_to_kg
    lbf_to_N = lbf_to_kg * g_SI
    pa_to_psi = (1 / (lbf_to_N)) * in2_to_m2
    bar_to_psi = pa_to_psi*100000

    # Energy
    BTU_to_J=1055.06
    kJ_to_BTU = 1000/BTU_to_J
    BTU_to_kJ = 1/kJ_to_BTU
    kJperkg_to_BTUperlb = kJ_to_BTU/kg_to_lbf
    m3perkg_to_ft3perlb = m3_to_ft3/kg_to_lbf

    @classmethod  # a classmethod can be used directly from a class without needing to instantiate an object
    def viscosityEnglishToSI(cls, mu, toSI=True):
        """
        Converts between lb*s/ft^2 and Pa*s
        :param mu: the viscosity in english units
        :param toSI:  True assumes english in, False assumes SI in
        :return: the viscosity in Pa*s if toSI=True, lb*s/ft^2 if toSI=False
        """
        # (lb*s)/ft^2*((3.3 ft/m)^2)*(1kg/2.2lb)*(9.81m/s^2)->(Pa*s)
        cf = (1 / cls.ft2_to_m2) * (cls.lbf_to_kg) * cls.g_SI
        return mu * cf if toSI else mu / cf

    @classmethod  # a classmethod can be used directly from a class without needing to instantiate an object
    def densityEnglishToSI(cls, rho, toSI=True):
        """
        Converts between lb/ft^3 and kg/m^3
        :param rho: specific weight or density
        :param toSI:  True assumes english in, False assumes SI in
        :return: density in SI or EN
        """
        # (lb/ft^3)*((3.3ft/m)^3)*(1kg/2.2lb) -> kg/m^3
        cf = cls.lbf_to_kg / cls.ft3_to_m3
        return rho * cf if toSI else rho / cf

    @classmethod  # a classmethod can be used directly from a class without needing to instantiate an object
    def head_to_pressure(cls, h, rho, SI=True):
        """
        Convert from height of column of fluid to pressure in consistent units
        :param h: head in height of fluid (in or m)
        :return: pressure in (psi or Pa)
        """
        if SI:  # p = rho*g*h = g*cf
            cf = rho * cls.g_SI / cls.gc_SI  # kg*m/m^3*s^2
            return h * cf
        else:  # p = rho*g*h = g*cf (h in in)
            cf = rho * cls.g_EN / cls.gc_EN * (1 / 12) ** 2  # (lbm*ft/ft^3*s^2)(lbf*s^2/lbm*ft)(ft^2/in^2)
            return h * cf
        # convert m of water to psi
        # (m)*(3.3*12in/m)*rho(kg/m^3)*(2.2lb/kg)*(1m/(3.3*12in))^3
        psi = p * cls.rho * 2.2 / ((3.3 * 12) ** 2)
        return psi

    @classmethod  # a classmethod can be used directly from a class without needing to instantiate an object
    def m_to_psi(cls, h, rho):
        """
        For converting from height of fluid to psi
        :param h: height of fluid in m
        :param rho: density of fluid in kg/m^3
        :return: pressure in psi
        """
        return cls.head_to_pressure(h, rho) * cls.pa_to_psi

    @classmethod  # a classmethod can be used directly from a class without needing to instantiate an object
    def psi_to_m(cls, p, rho):
        """
        For converting from psi to height of fluid.
        first convert psi to pa
        :param p: pressure in psi
        :param rho: density of fluid in kg/m^3
        :return: height of fluid in m
        """
        pa = p / cls.pa_to_psi
        h = pa / (rho * cls.g_SI)
        return h

    @classmethod # a classmethod can be used directly from a class without needing to instantiate an object
    def C_to_F(cls, T=0):
        return 9/5*(T)+32

    @classmethod # a classmethod can be used directly from a class without needing to instantiate an object
    def F_to_C(cls, T=0):
        return 5/9*(T-32)

    @classmethod # a classmethod can be used directly from a class without needing to instantiate an object
    def K_to_R(cls, T=0):
        return cls.C_to_F(T-273.15)+459.67

class SatPropsIsobar():
    def __init__(self, P):
        """
        Sets saturation properties for a given isobar
        :param P:  in kPa
        """
        tscol, pscol, hfcol, hgcol, sfcol, sgcol, vfcol, vgcol = np.loadtxt('sat_water_table.txt', skiprows=1, unpack=True)
        self.PSat=P/100
        self.TSat=float(griddata(pscol,tscol,self.PSat, method='cubic'))
        self.hf=float(griddata(pscol, hfcol,self.PSat,method = 'cubic'))
        self.hg=float(griddata(pscol,hgcol,self.PSat,method = 'cubic'))
        self.sf=float(griddata(pscol, sfcol,self.PSat,method = 'cubic'))
        self.sg=float(griddata(pscol,sgcol,self.PSat,method = 'cubic'))
        self.vf=float(griddata(pscol, vfcol,self.PSat,method = 'cubic'))
        self.vg=float(griddata(pscol,vgcol,self.PSat,method = 'cubic'))
        # for doing unit conversions
        self.UC=UnitConverter()
        self.getTextOutput()

    def getTextOutput(self, SI=True):
        """
        Sets the self.txtOut string for display.
        :param SI:
        :return: self.txtOut
        """
        if SI is False:
            P = self.PSat*self.UC.bar_to_psi
            PUnits="psi"
            T = self.UC.C_to_F(self.TSat)
            TUnits="F"
            hf=self.hf*self.UC.kJperkg_to_BTUperlb
            hg=self.hg*self.UC.kJperkg_to_BTUperlb
            HUnits="BTU/lb"
            sf=self.sf*self.UC.kJperkg_to_BTUperlb/self.UC.K_to_R(T=1)
            sg=self.sg*self.UC.kJperkg_to_BTUperlb/self.UC.K_to_R(T=1)
            SUnits="BTU/lb*R"
            vf=self.vf*self.UC.m3perkg_to_ft3perlb
            vg=self.vg*self.UC.m3perkg_to_ft3perlb
            VUnits="ft^3/lb"
        else:
            P = self.PSat
            PUnits = "bar"
            T = self.TSat
            TUnits = "C"
            hf = self.hf
            hg = self.hg
            HUnits = "kJ/kg"
            sf = self.sf
            sg = self.sg
            SUnits = "kJ/kg*K"
            vf = self.vf
            vg = self.vg
            VUnits = "m^3/kg"

        self.txtOut = "PSat = {:0.2f} {}, TSat = {:0.2f} {}".format(P, PUnits, T, TUnits)
        self.txtOut += "\nhf = {:0.2f} {}, hg = {:0.2f} {}".format(hf, HUnits,hg,HUnits)
        self.txtOut += "\nsf = {:0.2f} {}, sg = {:0.2f} {}".format(sf,SUnits,sg,SUnits)
        self.txtOut += "\nvf = {:0.4f} {}, vg = {:0.2f} {}".format(vf,VUnits,vg,VUnits)
        return self.txtOut

class Steam_SI:
    def __init__(self, P=None, T=None, x=None, v=None, h=None, s=None, name=None):
        """
        This is a general steam class for sub-critical (i.e., superheated and saturated) properties of steam.
        The user may specify any two properties to calculate all other properties of the steam.
        Note: we have 6 properties, but only can specify two of them.  Combinations=6!/(2!4!)=15

        I handle all cases in self.calc

        :param P: Pressure (kPa)
        :param T: Temperature (C)
        :param x: Quality
        :param v: Specific Volume (kg/m^3)
        :param h: Enthalpy (kJ/kg)
        :param s: Entropy (kJ/(kg*K))
        :param name:
        """
        self.P = P  # pressure - kPa
        self.T = T  # Temperature - degrees C
        self.x = x  # quality (a value between 0 and 1)
        self.v = v  # specific volume - m^3/kg
        self.h = h  # enthalpy - kJ/kg
        self.s = s  # entropy - kJ/(kg K)
        self.name = name # a useful identifier
        self.region = None # 'superheated' or 'saturated'
        self.RW=XSteam.specificGasConstant(self)
        self.calc()
        pass

    def calc(self):
        """
        In principle, there are 15 cases to handle any pair of variables.  I depend on user to specify proper set of variables
        :return: True if properties calculated
        """
        #region select case
        case=None
        if self.P is not None:  # pressure is specified
            if self.T is not None:
                case="PT"
            if self.x is not None:
                case="Px"
            if self.v is not None:
                case="Pv"
            if self.h is not None:
                case="Ph"
            if self.s is not None:
                case="Ps"
        if case is None and self.T is not None: # temperature is specified
            if self.x is not None:
                case="Tx"
            if self.v is not None:
                case="Tv"
            if self.h is not None:
                case="Th"
            if self.s is not None:
                case="Ts"
        if case is None and self.x is not None:  # quality is specified
            if self.v is not None:
                case ="xv"
            if self.h is not None:
                case="xh"
            if self.s is not None:
                case="xs"
        if case is None and self.v is not None:  # quality is specified
            if self.h is not None:
                case="vh"
            if self.s is not None:
                case="vs"
        if case is None and self.h is not None:  # enthalpy is specified
            if self.s is not None:
                case="hs"
        if case is None:
            return False
        # my 15 cases are:  PT, Px, Pv, Ph, Ps, Tx, Tv, Th, Ts, xv, xh, xs, vh, vs, hs
        #endregion

        st=XSteam(XSteam.UNIT_SYSTEM_MKS)  # m/kg/sec/°C/bar/W
        # Read the saturated table data
        tscol, pscol, hfcol, hgcol, sfcol, sgcol, vfcol, vgcol = np.loadtxt('sat_water_table.txt', skiprows=1, unpack=True)
        # Read the superheated water table
        tcol, hcol, scol, pcol = np.loadtxt('superheated_water_table.txt', skiprows=1, unpack=True)

        # ideal gas law V=RT/P-> v=(R/MW)*(T+273)/(P)
        # calculate a column of specific volume for the superheated water table using ideal gas
        vcol=[self.RW*(tcol[nI]+273)/pcol[nI] for nI in range(len(pcol))]

        #region calculate properties based on case
        # if P given (easy to check saturated vs. superheated)
        if case.__contains__("P"):
            # at the known pressure, get saturated liq and vap properties
            sfval = float(griddata(pscol, sfcol, self.P / 100, method='cubic'))
            sgval = float(griddata(pscol, sgcol, self.P / 100, method='cubic'))
            vfval = float(griddata(pscol, vfcol, self.P / 100, method='cubic'))
            vgval = float(griddata(pscol, vgcol, self.P / 100, method='cubic'))
            hfval = float(griddata(pscol, hfcol, self.P / 100, method='cubic'))
            hgval = float(griddata(pscol, hgcol, self.P / 100, method='cubic'))
            tsat = float(griddata(pscol, tscol, self.P / 100, method='cubic'))

            # this is a test of pyXSteam.  With cubic interpolation on tables, pretty good match
            sL=st.sL_p(self.P/100)
            sV=st.sV_p(self.P/100)
            vL=st.vL_p(self.P/100)
            vV=st.vV_p(self.P/100)
            hL=st.hL_p(self.P/100)
            hV=st.hV_p(self.P/100)
            t_sat=st.tsat_p(self.P/100)

            if case.__contains__("x"):  #quality given

                self.T=tsat
                self.s=self.x*(sgval-sfval)+sfval
                self.h=self.x*(hgval-hfval)+hfval
                self.v=self.x*(vgval-vfval)+vfval
                self.region="saturated"
                return True
            if case.__contains__("v"):  #find if saturated or superheated
                xval=(self.v-vfval)/(vgval-vfval)
                if xval<=1:  #saturated
                    self.region = "saturated"
                    self.x=xval
                    self.T = tsat
                    self.h = self.x * (hgval - hfval) + hfval
                    self.s = self.x * (sgval - sfval) + sfval
                    return True
                else:  #superheated
                    self.region = "superheated"
                    self.x=1.0
                    self.T = float(griddata((vcol, pcol), tcol, (self.v, self.P), method='cubic'))
                    self.h = float(griddata((vcol, pcol), hcol, (self.v, self.P), method='cubic'))
                    self.s = float(griddata((vcol, pcol), scol, (self.v, self.P), method='cubic'))
                    return True
            if case.__contains__("h"):  #find if saturated or superheated
                xval=(self.h-hfval)/(hgval-hfval)
                if xval<=1:  #saturated
                    self.region = "saturated"
                    self.x=xval
                    self.T = tsat
                    self.v = self.x * (vgval - vfval) + vfval
                    self.s = self.x * (sgval - sfval) + sfval
                    return True
                else:  #superheated
                    self.region = "superheated"
                    self.x=1.0
                    self.T = float(griddata((hcol, pcol), tcol, (self.h, self.P), method='cubic'))
                    self.v = self.igl_v()
                    self.s = float(griddata((hcol, pcol), scol, (self.h, self.P), method='cubic'))
                    return True
            if case.__contains__("s"):  # find if saturated or superheated
                xval = (self.s - sfval) / (sgval - sfval)
                if xval <= 1:  # saturated
                    self.region = "saturated"
                    self.x = xval
                    self.T = tsat
                    self.h = self.x * (hgval - hfval) + hfval
                    self.v = self.x * (vgval - vfval) + vfval
                    return True
                else:  # superheated
                    self.region = "superheated"
                    self.x=1.0
                    self.T = float(griddata((scol, pcol), tcol, (self.s, self.P), method='cubic'))
                    self.h = float(griddata((scol, pcol), hcol, (self.s, self.P), method='cubic'))
                    self.v = float(griddata((scol, pcol), vcol, (self.s, self.P), method='cubic'))
                    return True
            if case.__contains__("T"):  # find if satruated or superheated steam
                if self.T==tsat:  # generally, this is an indeterminate case, but assume saturated vapor
                    self.region = "saturated"
                    self.x=1
                    self.h=hgval
                    self.s=sgval
                    self.v=vgval
                    return True
                elif self.T>tsat:  #superheated
                    self.region = "superheated"
                    self.x=1.0
                    pbar=self.P/100

                    self.h=st.h_pt(self.P/100,self.T)
                    self.s=st.s_pt(self.P/100, self.T)
                    self.v = self.igl_v()
                    #self.h = float(griddata((tcol, pcol), hcol, (self.T, self.P), method='cubic'))
                    #self.s = float(griddata((tcol, pcol), scol, (self.T, self.P), method='cubic'))
                    return True
                else: # sub-cooled, so estimate properties
                    self.region = "saturated"
                    self.x=0
                    psat= float(griddata(tscol, pscol, self.T))*100
                    self.h=hfval+(self.P-psat)*vfval
                    self.s=sfval
                    self.v=vfval
                    return True
        # if T given (easy to check saturated vs. superheated)
        if case.__contains__("T"):
            # Using the known Temperature, interpolate on the saturation tables columns
            # at the known pressure
            sfval = float(griddata(tscol, sfcol, self.T, method='cubic'))
            sgval = float(griddata(tscol, sgcol, self.T, method='cubic'))
            vfval = float(griddata(tscol, vfcol, self.T, method='cubic'))
            vgval = float(griddata(tscol, vgcol, self.T, method='cubic'))
            hfval = float(griddata(tscol, hfcol, self.T, method='cubic'))
            hgval = float(griddata(tscol, hgcol, self.T, method='cubic'))
            psat = float(griddata(tscol, pscol, self.T, method='cubic'))
            if case.__contains__("x"):  # quality given
                self.region = "saturated"
                self.P = psat*100
                self.s = self.x * (sgval - sfval) + sfval
                self.h = self.x * (hgval - hfval) + hfval
                self.v = self.x * (vgval - vfval) + vfval
                return True
            if case.__contains__("v"):  # find if saturated or superheated
                xval = (self.v - vfval) / (vgval - vfval)
                if xval <= 1:  # saturated
                    self.region = "saturated"
                    self.x = xval
                    self.P = psat*100  #from bar to kPa
                    self.h = self.x * (hgval - hfval) + hfval
                    self.s = self.x * (sgval - sfval) + sfval
                    return True
                else:  # superheated
                    self.region = "superheated"
                    self.x=1.0
                    self.P = float(griddata((vcol, tcol), pcol, (self.v, self.T), method='cubic'))
                    self.h = float(griddata((vcol, tcol), hcol, (self.v, self.T), method='cubic'))
                    self.s = float(griddata((vcol, tcol), scol, (self.v, self.T), method='cubic'))
                    return True
            if case.__contains__("h"):  # find if saturated or superheated
                xval = (self.h - hfval) / (hgval - hfval)
                if xval <= 1:  # saturated
                    self.region = "saturated"
                    self.x = xval
                    self.P = psat*100
                    self.s = self.x * (sgval - sfval) + sfval
                    self.v = self.x * (vgval - vfval) + vfval
                    return True
                else:  # superheated
                    self.region = "superheated"
                    self.x=1.0
                    self.P = float(griddata((hcol, tcol), pcol, (self.h, self.T), method='cubic'))
                    self.s = float(griddata((hcol, tcol), scol, (self.h, self.T), method='cubic'))
                    self.v = self.igl_v()
                    return True
            if case.__contains__("s"):  # find if saturated or superheated
                xval = (self.s - sfval) / (sgval - sfval)
                if xval <= 1:  # saturated
                    self.region = "saturated"
                    self.x = xval
                    self.P = psat*100
                    self.h = self.x * (hgval - hfval) + hfval
                    self.v = self.x * (vgval - vfval) + vfval
                    return True
                else:  # superheated
                    self.region = "superheated"
                    self.x=1.0

                    self.P = float(griddata((scol, tcol), pcol, (self.s, self.T), method='cubic'))
                    self.h = float(griddata((scol, tcol), hcol, (self.s, self.T), method='cubic'))
                    self.v = self.igl_v()
                    return True
            return
        # if quality given (easy case) for saturated
        if case.__contains__("x"):
            self.region = "saturated"
            if case.__contains__("v"):  #find if saturated or superheated
                def findPsat_v(P):  #note: P in bar
                    vfval = float(griddata(pscol, vfcol, P, method='cubic'))
                    vgval = float(griddata(pscol, vgcol, P, method='cubic'))
                    vcalc=self.x*(vgval-vfval)+vfval
                    diff=self.v-vcalc
                    return diff
                self.P=100*fsolve(findPsat_v,[1])[0]  # find the pressure
                sfval = float(griddata(pscol, sfcol, self.P / 100, method='cubic'))
                sgval = float(griddata(pscol, sgcol, self.P / 100, method='cubic'))
                vfval = float(griddata(pscol, vfcol, self.P / 100, method='cubic'))
                vgval = float(griddata(pscol, vgcol, self.P / 100, method='cubic'))
                hfval = float(griddata(pscol, hfcol, self.P / 100, method='cubic'))
                hgval = float(griddata(pscol, hgcol, self.P / 100, method='cubic'))
                tsat = float(griddata(pscol, tscol, self.P / 100, method='cubic'))
                self.T = tsat
                self.s = self.x * (sgval - sfval) + sfval
                self.h = self.x * (hgval - hfval) + hfval
                return True
            if case.__contains__("h"):  #find if saturated or superheated
                def findPsat_h(P):  #note: P in bar
                    hfval = float(griddata(pscol, hfcol, P))
                    hgval = float(griddata(pscol, hgcol, P))
                    hcalc=self.x*(hgval-hfval)+hfval
                    return self.h-hcalc
                self.P=100*fsolve(findPsat_h,[1])[0]  # find the pressure
                sfval = float(griddata(pscol, sfcol, self.P / 100, method='cubic'))
                sgval = float(griddata(pscol, sgcol, self.P / 100, method='cubic'))
                vfval = float(griddata(pscol, vfcol, self.P / 100, method='cubic'))
                vgval = float(griddata(pscol, vgcol, self.P / 100, method='cubic'))
                hfval = float(griddata(pscol, hfcol, self.P / 100, method='cubic'))
                hgval = float(griddata(pscol, hgcol, self.P / 100, method='cubic'))
                tsat = float(griddata(pscol, tscol, self.P / 100, method='cubic'))
                self.T = tsat
                self.s = self.x * (sgval - sfval) + sfval
                self.v = self.x * (vgval - vfval) + vfval
                return True
            if case.__contains__("s"):  # find if saturated or superheated
                def findPsat_s(P):  #note: P in bar
                    sfval = float(griddata(pscol, sfcol, P, method='cubic'))
                    sgval = float(griddata(pscol, sgcol, P, method='cubic'))
                    stest=self.x*(sgval-sfval)+sfval
                    return self.s-stest
                self.P=100*fsolve(findPsat_s,[1])[0]  # find the pressure
                sfval = float(griddata(pscol, sfcol, self.P / 100, method='cubic'))
                sgval = float(griddata(pscol, sgcol, self.P / 100, method='cubic'))
                vfval = float(griddata(pscol, vfcol, self.P / 100, method='cubic'))
                vgval = float(griddata(pscol, vgcol, self.P / 100, method='cubic'))
                hfval = float(griddata(pscol, hfcol, self.P / 100, method='cubic'))
                hgval = float(griddata(pscol, hgcol, self.P / 100, method='cubic'))
                tsat = float(griddata(pscol, tscol, self.P / 100, method='cubic'))
                self.T = tsat
                self.h = self.x * (hgval - hfval) + hfval
                self.v = self.x * (vgval - vfval) + vfval
                return True
            return
        # if vh, vs, or hs (searching required to determine if saturated)
        if case.__contains__("v"):
            vcritical=float(griddata(pscol, vfcol, max(pscol), method='cubic'))
            twophase=False
            if self.v<vcritical:
                twophase=True
            else:  #might be two-phase
                tmax = float(griddata(vgcol, tscol, self.v, method='cubic'))  # temperature for x=1
                if case.__contains__("h"):
                    hgsat= float(griddata(tscol, hgcol, tmax, method='cubic'))
                    if self.h<=hgsat:  # means two phase
                        twophase=True
                else:
                    sgsat= float(griddata(tscol, sgcol, tmax, method='cubic'))
                    if self.s<=sgsat:  # means two phase
                        twophase=True
            if twophase:
                self.region = "saturated"
                if case.__contains__("h"):
                    tmax=max(tscol)
                    def findTSat_h(T):
                        hfsat=float(griddata(tscol, hfcol, T, method='cubic'))
                        hgsat=float(griddata(tscol, hgcol, T, method='cubic'))
                        vfsat=float(griddata(tscol, vfcol, T, method='cubic'))
                        vgsat=float(griddata(tscol, vgcol, T, method='cubic'))
                        xv=(self.v-vfsat)/(vgsat-vfsat)
                        xh=(self.h-hfsat)/(hgsat-hfsat)
                        return xh-xv if T<tmax else 100
                    self.T=fsolve(findTSat_h, 100)[0]
                    sfval = float(griddata(tscol, sfcol, self.T, method='cubic'))
                    sgval = float(griddata(tscol, sgcol, self.T, method='cubic'))
                    vfval = float(griddata(tscol, vfcol, self.T, method='cubic'))
                    vgval = float(griddata(tscol, vgcol, self.T, method='cubic'))
                    psat = float(griddata(tscol, pscol, self.T, method='cubic'))
                    self.x=(self.v-vfval)/(vgval-vfval)
                    self.P=psat*100
                    self.s=self.x*(sgval-sfval)+sfval
                    return True
                if case.__contains__("s"):
                    tmax=max(tscol)
                    def findTSat_s(T):
                        sfsat=float(griddata(tscol, sfcol, T, method='cubic'))
                        sgsat=float(griddata(tscol, sgcol, T, method='cubic'))
                        vfsat=float(griddata(tscol, vfcol, T, method='cubic'))
                        vgsat=float(griddata(tscol, vgcol, T, method='cubic'))
                        xv=(self.v-vfsat)/(vgsat-vfsat)
                        xs=(self.s-sfsat)/(sgsat-sfsat)
                        return xs-xv if T<tmax else 100
                    self.T=fsolve(findTSat_s, 100)[0]
                    vfval = float(griddata(tscol, vfcol, self.T, method='cubic'))
                    vgval = float(griddata(tscol, vgcol, self.T, method='cubic'))
                    hfval = float(griddata(tscol, hfcol, self.T, method='cubic'))
                    hgval = float(griddata(tscol, hgcol, self.T, method='cubic'))
                    psat = float(griddata(tscol, pscol, self.T, method='cubic'))
                    self.x=(self.v-vfval)/(vgval-vfval)
                    self.P=psat*100
                    self.h=self.x*(hgval-hfval)+hfval
                    return True
            else:  # means superheated
                self.region = "superheated"
                self.x = 1.0
                if case.__contains__("h"):
                    self.P = float(griddata((vcol, hcol), pcol, (self.v, self.h), method='cubic'))
                    self.T = float(griddata((vcol, hcol), tcol, (self.v, self.h), method='cubic'))
                    self.s = float(griddata((vcol, hcol), scol, (self.v, self.h), method='cubic'))
                    return True
                if case.__contains__("s"):
                    self.P = float(griddata((vcol, scol), pcol, (self.v, self.s), method='cubic'))
                    self.T = float(griddata((vcol, scol), tcol, (self.v, self.s), method='cubic'))
                    self.h = float(griddata((vcol, scol), hcol, (self.v, self.s), method='cubic'))
                    return True
                return True
        if case.__contains__("h"):  # this case has h & s given
            #see Mollier diagram.
            hV=float(griddata(sgcol, hgcol, self.s, method='cubic'))
            tmax=max(tscol)
            twophase=False
            if self.h<=hV:
                self.region="saturated"
                def findTSat_s(T):
                    sfsat=float(griddata(tscol, sfcol, T, method='cubic'))
                    sgsat=float(griddata(tscol, sgcol, T, method='cubic'))
                    hfsat=float(griddata(tscol, hfcol, T, method='cubic'))
                    hgsat=float(griddata(tscol, hgcol, T, method='cubic'))
                    xs=(self.s-sfsat)/(sgsat-sfsat)
                    xh=(self.h-hfsat)/(hgsat-hfsat)
                    return xs-xh if T<tmax else 100
                self.T=fsolve(findTSat_s, 100)[0]
                vfval = float(griddata(tscol, vfcol, self.T, method='cubic'))
                vgval = float(griddata(tscol, vgcol, self.T, method='cubic'))
                hfval = float(griddata(tscol, hfcol, self.T, method='cubic'))
                hgval = float(griddata(tscol, hgcol, self.T, method='cubic'))
                psat = float(griddata(tscol, pscol, self.T, method='cubic'))
                self.x=(self.h-hfval)/(hgval-hfval)
                self.P=psat*100
                self.v=self.x*(vgval-vfval)+vfval
                return True
            else:  # means superheated
                self.region="superheated"
                self.x=1.0
                self.P = float(griddata((hcol, scol), pcol, (self.h, self.s), method='cubic'))
                self.T = float(griddata((hcol, scol), tcol, (self.h, self.s), method='cubic'))
                self.v = self.igl_v()
                return True
        #endregion

    def igl_v(self):
        # ideal gas law V=RT/P-> v=(R/MW)*(T+273)/(P)
        # calculate a column of specific volume for the superheated water table using ideal gas
        return self.RW * (self.T + 273) / self.P

    def print(self):
        if self.name is not None:
            print('Name: {}'.format(self.name))
        if self.region is not None:
            print('Region: {}'.format(self.region))   
        if self.P is not None:
            print('p = {:.2f} kPa'.format(self.P))
        if self.T is not None:
            print('T = {:.1f} degrees C'.format(self.T))
        if self.h is not None:
            print('h = {:.2f} kJ/kg'.format(self.h))
        if self.s is not None:
            print('s = {:.4f} kJ/(kg K)'.format(self.s))
        if self.v is not None:
            print('v = {:.6f} m^3/kg'.format(self.v))
        if self.x is not None:
            print('x = {:.4f}'.format(self.x))
        print('')

def main():
    inlet=Steam_SI(P=7350, name='Turbine Inlet')
    inlet.x =0.9  # 90 percent quality  - accessed DIRECTLY  - not through a method
    inlet.calc()
    inlet.print()
    h1 = inlet.h; s1 = inlet.s
    print(h1,s1,'\n')

    outlet=Steam_SI(P=100, s=inlet.s, name='Turbine Exit')
           #notice -  s=inlet.s
    outlet.print()


    another=Steam_SI(P=8575, h=2050, name='State 3')
    another.print()

    yetanother=Steam_SI(P=8900, h=41250, name='State 4')
    yetanother.print()

    g1=Steam_SI(P=800, x=1.0, name='Gap1')
    g1.print()
    g2=Steam_SI(P=g1.P, s=g1.s * 1.0012, name='Gap2')
    g2.print()
    g2=Steam_SI(P=g1.P, s=6.6699, name='Gap3')
    g2.print()




    #uncommenting the next two lines will cause a ValueError error
    #final1 = State(7000, T=250, name='State 5')
    #final1.print()

if __name__ == "__main__":
   main()

