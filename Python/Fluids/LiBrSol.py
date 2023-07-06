import math
import numpy
import pandas  
import CoolProp.CoolProp as CP
from scipy.optimize import newton, root_scalar 
DFcal_LiBrSol = pandas.read_csv("Fluids/LiBrSol_cal.csv", sep=";")

def  Calc_cp_from_T_X_LiBrSol_Patek(T,X_mol_LiBr):
    """
    Calc_cp_from_T_X_LiBrSol_Patek
    Uses coefficients and formula obtained from Patek 2006
    Input:
        -   Temperature of Solution T                                   [K]
        -   Molar Concentration X of LiBr in Solution                   [-]
    Output:
        -   Molar heat capacity                                         [J/mol/K]
    """
    ## Constants
    cp_t = 76.0226             #[J/molK]
    T_c = 647.096              #[K]
    T_0 = 221                  #[K]
    RMS = 0.98 * 10**-3        #[-]
    T_t = 273.16               #[K]
    # Table 6
    Koef_a = [-1.42094  *   10**1,
            4.04943   *   10**1,
            1.11135   *   10**2,
            2.29980   *   10**2,
            1.34526   *   10**3,
            -1.41010  *   10**-2,
            1.24977   *   10**-2,
            -6.83209  *   10**-4]
    Koef_t = [0,0,0,0,0,2,3,4]
    Koef_n = [0,0,1,2,3,0,3,2]
    Koef_m = [2,3,3,3,3,2,1,1]
    # Table 13
    Koef_beta = [0,2,3,6,34]
    Koef_gamma = [0,2,3,5,0]
    Koef_alpha = [1.38801,
                -2.95318,
                3.18721,
                -0.645473,
                9.18946 * 10**5]
    ## Calculation
    # Calculation of cp_sat
    sum = 0
    for i in range(0,5,1):
        sum = sum + Koef_alpha[i]*(1-(T/T_c))**Koef_beta[i]*(T/T_t)**Koef_gamma[i]
    cp_sat = cp_t*sum
    # Calculation of cp
    cp = (1-X_mol_LiBr)*cp_sat+\
        cp_t*(Koef_a[0]*X_mol_LiBr**Koef_m[0]*(0.4 - X_mol_LiBr)**Koef_n[0]*(T_c/(T-T_0))**Koef_t[0]+\
        Koef_a[1]*X_mol_LiBr**Koef_m[1]*(0.4 - X_mol_LiBr)**Koef_n[1]*(T_c/(T-T_0))**Koef_t[1]+\
        Koef_a[2]*X_mol_LiBr**Koef_m[2]*(0.4 - X_mol_LiBr)**Koef_n[2]*(T_c/(T-T_0))**Koef_t[2]+\
        Koef_a[3]*X_mol_LiBr**Koef_m[3]*(0.4 - X_mol_LiBr)**Koef_n[3]*(T_c/(T-T_0))**Koef_t[3]+\
        Koef_a[4]*X_mol_LiBr**Koef_m[4]*(0.4 - X_mol_LiBr)**Koef_n[4]*(T_c/(T-T_0))**Koef_t[4]+\
        Koef_a[5]*X_mol_LiBr**Koef_m[5]*(0.4 - X_mol_LiBr)**Koef_n[5]*(T_c/(T-T_0))**Koef_t[5]+\
        Koef_a[6]*X_mol_LiBr**Koef_m[6]*(0.4 - X_mol_LiBr)**Koef_n[6]*(T_c/(T-T_0))**Koef_t[6]+\
        Koef_a[7]*X_mol_LiBr**Koef_m[7]*(0.4 - X_mol_LiBr)**Koef_n[7]*(T_c/(T-T_0))**Koef_t[7])
    return  cp


def Calc_h_from_T_X_LiBrSol_Patek(T,X_mol_LiBr):
    """
    Calc_h_from_T_X_LiBrSol_Patek
    Uses coefficients and formula obtained from Patek 2006
    Input:
        -   Temperature of Solution T                                   [K]
        -   Molar Concentration X of LiBr in Solution                   [-]
    Output:
        -   Molar enthalpy h of solution                                [J/mol]
    """
    ## Constants
    T_c = 647.096               #[K]
    h_c = 37548.5               #[J/mol]
    T_0 = 221                   #[K]
    # Table 7
    Koef_a=Koef_a = [  2.27431     *   10**0,
                -7.99511    *   10**0,
                3.85239     *   10**2,
                -1.63940    *   10**4,
                -4.22562    *   10**2,          
                1.13314     *   10**-1,
                -8.33474    *   10**0,
                -1.73833    *   10**4,
                6.49763     *   10**0,
                3.24552     *   10**3,
                -1.34643    *   10**4,
                3.99322     *   10**4,
                -2.58877    *   10**5,
                -1.93046    *   10**-3,
                2.80616     *   10**0,
                -4.04479    *   10**1,
                1.45342     *   10**2,
                -2.74873    *   10**0,
                -4.49743    *   10**2,
                -1.21794    *   10**1,
                -5.83739    *   10**-3,
                2.33910     *   10**-1,
                3.41888     *   10**-1,
                8.85259     *   10**0,
                -1.78731    *   10**1,
                7.35179     *   10**-2,
                -1.79430    *   10**-4,
                1.84261     *   10**-3,
                -6.24282    *   10**-3,
                6.84765     *   10**-3]
    Koef_t = [0,0,0,0,0,1,1,1,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5]
    Koef_n = [0,1,6,6,2,0,0,4,0,4,5,5,6,0,3,5,7,0,3,1,0,4,2,6,7,0,0,1,2,3]
    Koef_m = [1,1,2,3,6,1,3,5,4,5,5,6,6,1,2,2,2,5,6,7,1,1,2,2,2,3,1,1,1,1]
    # Table 14
    Koef_beta = [1/3,2/3,5/6,21/6]
    Koef_alpha = [  -4.37196    *   10**-1,
                    3.03440     *   10**-1,
                    -1.29582    *   10**0,
                    -1.76410    *   10**-1]
    ## Calculation
    #  Calculation of h_start
    sum = 0
    for i in range(4):
        sum = sum + Koef_alpha[i]*(1-(T/T_c))**Koef_beta[i]
    h_sat = h_c * (1 + sum)
    # Calculation of h
    factors = numpy.zeros(30)   #numpy.zeros((30),)
    a=0
    b=0
    c=0
    d=0
    e=0
    f=0
    for i in range(30):
        factors[i] =  Koef_a[i]*X_mol_LiBr**Koef_m[i]*(0.4-X_mol_LiBr)**Koef_n[i]
        if Koef_t[i] == 0:
            f = f + factors[i]
        elif Koef_t[i] == 1:
            e = e + factors[i]
        elif Koef_t[i] == 2:
            d = d + factors[i]
        elif Koef_t[i] == 3:
            c = c + factors[i]
        elif Koef_t[i] == 4:
            b = b + factors[i]
        elif Koef_t[i] == 5:
            a = a + factors[i]
    h = (1-X_mol_LiBr)*h_sat + h_c*(a*(T_c/(T-T_0))**5 + b*(T_c/(T-T_0))**4 + c*(T_c/(T-T_0))**3 + d*(T_c/(T-T_0))**2 + e*(T_c/(T-T_0))**1 + f)
    return h


def Calc_p_from_T_X_LiBrSol_Patek(T,X_mol_LiBr):
    """
    Calc_p_from_T_X_LiBrSol_Patek
    Uses coefficients and formula obtained from Patek 2006
    Input:
        -   Temperature of Solution T                                   [K]
        -   Molar Concentration X of LiBr in Solution                   [-]
    Output:
        -   Pressure of solution at saturation                          [Pa]
    """
    ## Constants
    T_c = 647.096              #[K]
    p_c = 22.064e6             #[Pa]
    # Table 4
    Koef_a = [  -2.41303    *   10**2,
                1.91750     *   10**7,
                -1.75521    *   10**8,
                3.25430     *   10**7,
                3.92571     *   10**2,          
                -2.12626    *   10**3,
                1.85127     *   10**8,
                1.91216     *   10**3]
    Koef_t = [0,0,0,0,1,1,1,1]
    Koef_n = [0,5,6,3,0,2,6,0]
    Koef_m = [3,4,4,8,1,1,4,6]
    # Table 11
    Koef_beta = [1.0,1.5,3.0,3.5,4.0,7.5]
    Koef_alpha = [  -7.85951783,
                    1.84408259,
                    -11.7866497,
                    22.6807411,
                    -15.9618719,
                    1.80122502]
    ## Calculation
    # Calculation of Teta
    sum = 0
    for i in range(8):
        sum = sum +  Koef_a[i]*X_mol_LiBr**Koef_m[i]*(0.4-X_mol_LiBr)**Koef_n[i]*(T/T_c)**Koef_t[i]
    Teta = T - sum
    # Calculation of vapor pressure p_vap
    sum = 0
    for i in range(6):
        sum = sum + Koef_alpha[i]*(1-(Teta/T_c))**Koef_beta[i]
    p_vap = p_c*math.exp((T_c/Teta)*sum)
    return p_vap


def Calc_rho_from_T_X_LiBrSol_Patek(T,X_mol_LiBr):
    """
    Calc_rho_from_T_X_LiBrSol_Patek
    Uses coefficients and formula obtained from Patek 2006
    Input:
        -   Temperature of Solution T                                   [K]
        -   Molar Concentration X of LiBr in Solution                   [-]
    Output:
        -   Molar density rho                                           [mol/m**3]
    """
    ## Constants
    T_c = 647.096             #[K]
    rho_c = 17.873             #[mol/m**3]
    RMS = 0.44 * 10**-3         #[-]
    # Table 5
    Koef_a = [1.746,4.709]
    Koef_t = [0,6]
    Koef_m = [1,1]
    # Table 12
    Koef_beta = [1/3,2/3,5/3,16/3,43/3,110/3]
    Koef_alpha = [  1.99274064,
                    1.09965342,
                    -0.510839303,
                    -1.75493479,
                    -45.5170352,
                    -6.7469445  *   10**5]
    ## Calculation
    # Calculation of rho_sat
    sum = 0
    for i in range(6):
        sum = sum + Koef_alpha[i]*(1-(T/T_c))**Koef_beta[i]
    rho_sat = rho_c * (1 + sum)
    # Calculation of rho
    rho = (1-X_mol_LiBr)*rho_sat+\
    rho_c*(Koef_a[0]*X_mol_LiBr**Koef_m[0]*(T/T_c)**Koef_t[0]+\
    Koef_a[1]*X_mol_LiBr**Koef_m[1]*(T/T_c)**Koef_t[1])
    rho = rho*1000
    return rho


def Calc_s_from_T_X_LiBrSol_Patek(T,X_mol_LiBr):
    """
    Calc_s_from_T_X_LiBrSol_Patek
    Uses coefficients and formula obtained from Patek 2006
    Input:
        -   Temperature of Solution T                                   [K]
        -   Molar concentration of LiBr in Solution                     [-]
    Output:
        -   Molar entropy of solution                                   [J/mol/K]
    """
    ## Constants
    T_c = 647.096              #[K]
    s_c = 79.3933              #[J/molK]
    T_0 = 221                  #[K]
    # Table 8
    Koef_a = [  1.53091     *   10**0,
                -4.52564    *   10**0,
                6.98302     *   10**2,
                -2.1666     *   10**4,
                -1.47533    *   10**3,
                8.47012     *   10**-2,
                -6.59523    *   10**0,
                -2.95331    *   10**4,
                9.56314     *   10**-3,
                -1.88679    *   10**-1,
                9.31752     *   10**0,
                5.78104     *   10**0,
                1.38931     *   10**4,
                -1.71762    *   10**4,
                4.15108     *   10**2,
                -5.55647    *   10**4,
                -4.23409    *   10**-3,
                3.05242     *   10**1,
                -1.67620    *   10**0,
                1.48283     *   10**1,
                3.03055     *   10**-3,
                -4.01810    *   10**-2,
                1.49252     *   10**-1,
                2.59240     *   10**0,
                -1.77421    *   10**-1,
                -6.99650    *   10**-5,
                6.05007     *   10**-4,
                -1.65228    *   10**-3,
                1.22966     *   10**-3]
    Koef_t = [0,0,0,0,0,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,4,4,5,5,5,5]
    Koef_n = [0,1,6,6,2,0,0,4,0,0,4,0,4,5,2,5,0,4,0,1,0,2,4,7,1,0,1,2,3]
    Koef_m = [1,1,2,3,6,1,3,5,1,2,2,4,5,5,6,6,1,3,5,7,1,1,1,2,3,1,1,1,1]
    # Table 15
    Koef_beta = [1/3,1,8/3,8]
    Koef_alpha = [  -3.34112    *   10**-1,
                    -8.47987    *   10**-1,
                    -9.11980    *   10**-1,
                    -1.64046    *   10**0]
    ## Calculation
    # Calculation of s_sat
    sum = 0
    for i in range(4):
        sum = sum + Koef_alpha[i]*(1-(T/T_c))**Koef_beta[i]
    s_sat = s_c * (1 + sum)
    # Calculation of s
    factors = numpy.zeros((29,))
    a = 0
    b = 0
    c = 0
    d = 0
    e = 0
    f = 0
    for i in range(29):
        factors[i] =  Koef_a[i]*X_mol_LiBr**Koef_m[i]*(0.4-X_mol_LiBr)**Koef_n[i]
        if Koef_t[i] == 0:
            f = f + factors[i]
        elif Koef_t[i] == 1:
            e = e + factors[i]
        elif Koef_t[i] == 2:
            d = d + factors[i]
        elif Koef_t[i] == 3:
            c = c + factors[i]
        elif Koef_t[i] == 4:
            b = b + factors[i]
        elif Koef_t[i] == 5:
            a = a + factors[i]
    s = (1-X_mol_LiBr)*s_sat + s_c*(a*(T_c/(T-T_0))**5 + b*(T_c/(T-T_0))**4 + c*(T_c/(T-T_0))**3 + d*(T_c/(T-T_0))**2 + e*(T_c/(T-T_0))**1 + f)
    return s


def  Calc_T_from_h_X_LiBrSol_Patek(h,X_mol_LiBr):
    """
    Calc_T_from_h_X_LiBrSol_Patek
    Uses coefficients and formula obtained from Patek 2006
    Input:
        -   Molar enthalpy of solution                                      [J/mol/K]
        -   Molar concentration X of LiBr in Solution                       [-]
    Output:
        -   Temperature of solution                                         [K]
    """
    #Constants
    T_c = 647.096              #[K]
    h_c = 37548.5              #[J/mol]
    T_0 = 221                  #[K]
    #Table 7
    Koef_a =[2.27431 * 10**0,
                -7.99511    *   10**0,
                3.85239     *   10**2,
                -1.63940    *   10**4,
                -4.22562    *   10**2,          
                1.13314     *   10**-1,
                -8.33474    *   10**0,
                -1.73833    *   10**4,
                6.49763     *   10**0,
                3.24552     *   10**3,
                -1.34643    *   10**4,
                3.99322     *   10**4,
                -2.58877    *   10**5,
                -1.93046    *   10**-3,
                2.80616     *   10**0,
                -4.04479    *   10**1,
                1.45342     *   10**2,
                -2.74873    *   10**0,
                -4.49743    *   10**2,
                -1.21794    *   10**1,
                -5.83739    *   10**-3,
                2.33910     *   10**-1,
                3.41888     *   10**-1,
                8.85259     *   10**0,
                -1.78731    *   10**1,
                7.35179     *   10**-2,
                -1.79430    *   10**-4,
                1.84261     *   10**-3,
                -6.24282    *   10**-3,
                6.84765     *   10**-3]
    Koef_t = [0,0,0,0,0,1,1,1,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5]
    Koef_n = [0,1,6,6,2,0,0,4,0,4,5,5,6,0,3,5,7,0,3,1,0,4,2,6,7,0,0,1,2,3]
    Koef_m = [1,1,2,3,6,1,3,5,4,5,5,6,6,1,2,2,2,5,6,7,1,1,2,2,2,3,1,1,1,1]
    # Table 14
    Koef_beta = [1/3,2/3,5/6,21/6]
    Koef_alpha = [  -4.37196    *   10**-1,
                        3.03440     *   10**-1,
                        -1.29582    *   10**0,
                        -1.76410    *   10**-1]
    ## Calculation
    fun = lambda T :(1-X_mol_LiBr)*h_c*(1+\
            Koef_alpha(0)*(1-(T/T_c))**Koef_beta(0)+\
            Koef_alpha(1)*(1-(T/T_c))**Koef_beta(1)+\
            Koef_alpha(2)*(1-(T/T_c))**Koef_beta(2)+\
            Koef_alpha(3)*(1-(T/T_c))**Koef_beta(3))+\
            h_c*\
            (Koef_a(0)*X_mol_LiBr**Koef_m(0)*(0.4-X_mol_LiBr)**Koef_n(0)*(T_c/(T-T_0))**Koef_t(0)+\
            Koef_a(1)*X_mol_LiBr**Koef_m(1)*(0.4-X_mol_LiBr)**Koef_n(1)*(T_c/(T-T_0))**Koef_t(1)+\
            Koef_a(2)*X_mol_LiBr**Koef_m(2)*(0.4-X_mol_LiBr)**Koef_n(2)*(T_c/(T-T_0))**Koef_t(2)+\
            Koef_a(3)*X_mol_LiBr**Koef_m(3)*(0.4-X_mol_LiBr)**Koef_n(3)*(T_c/(T-T_0))**Koef_t(3)+\
            Koef_a(4)*X_mol_LiBr**Koef_m(4)*(0.4-X_mol_LiBr)**Koef_n(4)*(T_c/(T-T_0))**Koef_t(4)+\
            Koef_a(5)*X_mol_LiBr**Koef_m(5)*(0.4-X_mol_LiBr)**Koef_n(5)*(T_c/(T-T_0))**Koef_t(5)+\
            Koef_a(6)*X_mol_LiBr**Koef_m(6)*(0.4-X_mol_LiBr)**Koef_n(6)*(T_c/(T-T_0))**Koef_t(6)+\
            Koef_a(7)*X_mol_LiBr**Koef_m(7)*(0.4-X_mol_LiBr)**Koef_n(7)*(T_c/(T-T_0))**Koef_t(7)+\
            Koef_a(8)*X_mol_LiBr**Koef_m(8)*(0.4-X_mol_LiBr)**Koef_n(8)*(T_c/(T-T_0))**Koef_t(8)+\
            Koef_a(9)*X_mol_LiBr**Koef_m(9)*(0.4-X_mol_LiBr)**Koef_n(9)*(T_c/(T-T_0))**Koef_t(9)+\
            Koef_a(10)*X_mol_LiBr**Koef_m(10)*(0.4-X_mol_LiBr)**Koef_n(10)*(T_c/(T-T_0))**Koef_t(10)+\
            Koef_a(11)*X_mol_LiBr**Koef_m(11)*(0.4-X_mol_LiBr)**Koef_n(11)*(T_c/(T-T_0))**Koef_t(11)+\
            Koef_a(12)*X_mol_LiBr**Koef_m(12)*(0.4-X_mol_LiBr)**Koef_n(12)*(T_c/(T-T_0))**Koef_t(12)+\
            Koef_a(13)*X_mol_LiBr**Koef_m(13)*(0.4-X_mol_LiBr)**Koef_n(13)*(T_c/(T-T_0))**Koef_t(13)+\
            Koef_a(14)*X_mol_LiBr**Koef_m(14)*(0.4-X_mol_LiBr)**Koef_n(14)*(T_c/(T-T_0))**Koef_t(14)+\
            Koef_a(15)*X_mol_LiBr**Koef_m(15)*(0.4-X_mol_LiBr)**Koef_n(15)*(T_c/(T-T_0))**Koef_t(15)+\
            Koef_a(16)*X_mol_LiBr**Koef_m(16)*(0.4-X_mol_LiBr)**Koef_n(16)*(T_c/(T-T_0))**Koef_t(16)+\
            Koef_a(17)*X_mol_LiBr**Koef_m(17)*(0.4-X_mol_LiBr)**Koef_n(17)*(T_c/(T-T_0))**Koef_t(17)+\
            Koef_a(18)*X_mol_LiBr**Koef_m(18)*(0.4-X_mol_LiBr)**Koef_n(18)*(T_c/(T-T_0))**Koef_t(18)+\
            Koef_a(19)*X_mol_LiBr**Koef_m(19)*(0.4-X_mol_LiBr)**Koef_n(19)*(T_c/(T-T_0))**Koef_t(19)+\
            Koef_a(20)*X_mol_LiBr**Koef_m(20)*(0.4-X_mol_LiBr)**Koef_n(20)*(T_c/(T-T_0))**Koef_t(20)+\
            Koef_a(21)*X_mol_LiBr**Koef_m(21)*(0.4-X_mol_LiBr)**Koef_n(21)*(T_c/(T-T_0))**Koef_t(21)+\
            Koef_a(22)*X_mol_LiBr**Koef_m(22)*(0.4-X_mol_LiBr)**Koef_n(22)*(T_c/(T-T_0))**Koef_t(22)+\
            Koef_a(23)*X_mol_LiBr**Koef_m(23)*(0.4-X_mol_LiBr)**Koef_n(23)*(T_c/(T-T_0))**Koef_t(23)+\
            Koef_a(24)*X_mol_LiBr**Koef_m(24)*(0.4-X_mol_LiBr)**Koef_n(24)*(T_c/(T-T_0))**Koef_t(24)+\
            Koef_a(25)*X_mol_LiBr**Koef_m(25)*(0.4-X_mol_LiBr)**Koef_n(25)*(T_c/(T-T_0))**Koef_t(25)+\
            Koef_a(26)*X_mol_LiBr**Koef_m(26)*(0.4-X_mol_LiBr)**Koef_n(26)*(T_c/(T-T_0))**Koef_t(26)+\
            Koef_a(27)*X_mol_LiBr**Koef_m(27)*(0.4-X_mol_LiBr)**Koef_n(27)*(T_c/(T-T_0))**Koef_t(27)+\
            Koef_a(28)*X_mol_LiBr**Koef_m(28)*(0.4-X_mol_LiBr)**Koef_n(28)*(T_c/(T-T_0))**Koef_t(28)+\
            Koef_a(29)*X_mol_LiBr**Koef_m(29)*(0.4-X_mol_LiBr)**Koef_n(29)*(T_c/(T-T_0))**Koef_t(29))-h
    T = newton(fun,x0=300)
    return  T 


def  Calc_T_from_p_X_satLiBrSol_Patek(p,X):
    """
    Calc_T_from_p_X_satLiBrSol_Patek
    Uses coefficients and formula obtained from Patek 2006
    Input:
        -   pressure of LiBr solution                                  [Pa]
        -   concentration of  solution                                 [-]
    Output:
        -   Saturation temperature of LiBr solution                    [T]
    """
    ## Constants
    T_c = 647.096              #[K]
    p_c = 22.064 * 10**6        #[Pa]
    # Table 4
    Koef_a = [  -2.41303    *   10**2,
                1.91750     *   10**7,
                -1.75521    *   10**8,
                3.25430     *   10**7,
                3.92571     *   10**2,
                -2.12626    *   10**3,
                1.85127     *   10**8,
                1.91216     *   10**3]
    Koef_t = [0,0,0,0,1,1,1,1]
    Koef_n = [0,5,6,3,0,2,6,0]
    Koef_m = [3,4,4,8,1,1,4,6]
    # Table 11
    Koef_beta = [1.0,1.5,3.0,3.5,4.0,7.5]
    Koef_alpha = [  -7.85951783,
                    1.84408259,
                    -11.7866497,
                    22.6807411,
                    -15.9618719,
                    1.80122502]
    ## Calculation
    # Calculate T_sat
    funT = lambda T:p_c*math.exp((T_c/T)*\
        (Koef_alpha[0]*(1-(T/T_c))**Koef_beta[0]+\
        Koef_alpha[1]*(1-(T/T_c))**Koef_beta[1]+\
        Koef_alpha[2]*(1-(T/T_c))**Koef_beta[2]+\
        Koef_alpha[3]*(1-(T/T_c))**Koef_beta[3]+\
        Koef_alpha[4]*(1-(T/T_c))**Koef_beta[4]+\
        Koef_alpha[5]*(1-(T/T_c))**Koef_beta[5]))-p
    T_sat = newton(funT,x0=300)
    funT = lambda T: T-(Koef_a[0]*X**Koef_m[0]*(0.4-X)**Koef_n[0]*(T/T_c)**Koef_t[0]\
        +Koef_a[1]*X**Koef_m[1]*(0.4-X)**Koef_n[1]*(T/T_c)**Koef_t[1]\
        +Koef_a[2]*X**Koef_m[2]*(0.4-X)**Koef_n[2]*(T/T_c)**Koef_t[2]\
        +Koef_a[3]*X**Koef_m[3]*(0.4-X)**Koef_n[3]*(T/T_c)**Koef_t[3]\
        +Koef_a[4]*X**Koef_m[4]*(0.4-X)**Koef_n[4]*(T/T_c)**Koef_t[4]\
        +Koef_a[5]*X**Koef_m[5]*(0.4-X)**Koef_n[5]*(T/T_c)**Koef_t[5]\
        +Koef_a[6]*X**Koef_m[6]*(0.4-X)**Koef_n[6]*(T/T_c)**Koef_t[6]\
        +Koef_a[7]*X**Koef_m[7]*(0.4-X)**Koef_n[7]*(T/T_c)**Koef_t[7])-T_sat
    T = newton(funT,x0=300)
    #T=root_scalar(funT, bracket=[0,800])
    return  T 


def  Calc_T_from_s_X_LiBrSol_Patek(s,X_mol_LiBr):
    """
    Calc_T_from_s_X_LiBrSol_Patek
    Uses coefficients and formula obtained from Patek 2006
    Input:
        -   Molar enthalpy of solution                                  [J/mol]
        -   Molar Concentration X of LiBr in Solution                   [-]
    Output:
        -   Temperature of solution                                     [K]
    """
    ## Constants
    T_c = 647.096              #[K]
    s_c = 79.3933              #[J/molK]
    T_0 = 221                  #[K]
    # Table 8
    Koef_a = [  1.53091     *   10**0,
                -4.52564    *   10**0,
                6.98302     *   10**2,
                -2.1666     *   10**4,
                -1.47533    *   10**3,
                8.47012     *   10**-2,
                -6.59523    *   10**0,
                -2.95331    *   10**4,
                9.56314     *   10**-3,
                -1.88679    *   10**-1,
                9.31752     *   10**0,
                5.78104     *   10**0,
                1.38931     *   10**4,
                -1.71762    *   10**4,
                4.15108     *   10**2,
                -5.55647    *   10**4,
                -4.23409    *   10**-3,
                3.05242     *   10**1,
                -1.67620    *   10**0,
                1.48283     *   10**1,
                3.03055     *   10**-3,
                -4.01810    *   10**-2,
                1.49252     *   10**-1,
                2.59240     *   10**0,
                -1.77421    *   10**-1,
                -6.99650    *   10**-5,
                6.05007     *   10**-4,
                -1.65228    *   10**-3,
                1.22966     *   10**-3]
    Koef_t = [0,0,0,0,0,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,4,4,5,5,5,5]
    Koef_n = [0,1,6,6,2,0,0,4,0,0,4,0,4,5,2,5,0,4,0,1,0,2,4,7,1,0,1,2,3]
    Koef_m = [1,1,2,3,6,1,3,5,1,2,2,4,5,5,6,6,1,3,5,7,1,1,1,2,3,1,1,1,1]
    # Table 15
    Koef_beta = [1/3,1,8/3,8]
    Koef_alpha = [  -3.34112    *   10**-1,
                    -8.47987    *   10**-1,
                    -9.11980    *   10**-1,
                    -1.64046    *   10**0]
    ## Calculation
    fun = lambda T : (1-X_mol_LiBr)*s_c*(1+\
        Koef_alpha[0]*(1-(T/T_c))**Koef_beta[0]+\
        Koef_alpha[1]*(1-(T/T_c))**Koef_beta[1]+\
        Koef_alpha[2]*(1-(T/T_c))**Koef_beta[2]+\
        Koef_alpha[3]*(1-(T/T_c))**Koef_beta[3])+\
        s_c*\
        (Koef_a[0]*X_mol_LiBr**Koef_m[0]*(0.4-X_mol_LiBr)**Koef_n[0]*(T_c/(T-T_0))**Koef_t[0]+\
        Koef_a[1]*X_mol_LiBr**Koef_m[1]*(0.4-X_mol_LiBr)**Koef_n[1]*(T_c/(T-T_0))**Koef_t[1]+\
        Koef_a[2]*X_mol_LiBr**Koef_m[2]*(0.4-X_mol_LiBr)**Koef_n[2]*(T_c/(T-T_0))**Koef_t[2]+\
        Koef_a[3]*X_mol_LiBr**Koef_m[3]*(0.4-X_mol_LiBr)**Koef_n[3]*(T_c/(T-T_0))**Koef_t[3]+\
        Koef_a[4]*X_mol_LiBr**Koef_m[4]*(0.4-X_mol_LiBr)**Koef_n[4]*(T_c/(T-T_0))**Koef_t[4]+\
        Koef_a[5]*X_mol_LiBr**Koef_m[5]*(0.4-X_mol_LiBr)**Koef_n[5]*(T_c/(T-T_0))**Koef_t[5]+\
        Koef_a[6]*X_mol_LiBr**Koef_m[6]*(0.4-X_mol_LiBr)**Koef_n[6]*(T_c/(T-T_0))**Koef_t[6]+\
        Koef_a[7]*X_mol_LiBr**Koef_m[7]*(0.4-X_mol_LiBr)**Koef_n[7]*(T_c/(T-T_0))**Koef_t[7]+\
        Koef_a[8]*X_mol_LiBr**Koef_m[8]*(0.4-X_mol_LiBr)**Koef_n[8]*(T_c/(T-T_0))**Koef_t[8]+\
        Koef_a[9]*X_mol_LiBr**Koef_m[9]*(0.4-X_mol_LiBr)**Koef_n[9]*(T_c/(T-T_0))**Koef_t[9]+\
        Koef_a[10]*X_mol_LiBr**Koef_m[10]*(0.4-X_mol_LiBr)**Koef_n[10]*(T_c/(T-T_0))**Koef_t[10]+\
        Koef_a[11]*X_mol_LiBr**Koef_m[11]*(0.4-X_mol_LiBr)**Koef_n[11]*(T_c/(T-T_0))**Koef_t[11]+\
        Koef_a[12]*X_mol_LiBr**Koef_m[12]*(0.4-X_mol_LiBr)**Koef_n[12]*(T_c/(T-T_0))**Koef_t[12]+\
        Koef_a[13]*X_mol_LiBr**Koef_m[13]*(0.4-X_mol_LiBr)**Koef_n[13]*(T_c/(T-T_0))**Koef_t[13]+\
        Koef_a[14]*X_mol_LiBr**Koef_m[14]*(0.4-X_mol_LiBr)**Koef_n[14]*(T_c/(T-T_0))**Koef_t[14]+\
        Koef_a[15]*X_mol_LiBr**Koef_m[15]*(0.4-X_mol_LiBr)**Koef_n[15]*(T_c/(T-T_0))**Koef_t[15]+\
        Koef_a[16]*X_mol_LiBr**Koef_m[16]*(0.4-X_mol_LiBr)**Koef_n[16]*(T_c/(T-T_0))**Koef_t[16]+\
        Koef_a[17]*X_mol_LiBr**Koef_m[17]*(0.4-X_mol_LiBr)**Koef_n[17]*(T_c/(T-T_0))**Koef_t[17]+\
        Koef_a[18]*X_mol_LiBr**Koef_m[18]*(0.4-X_mol_LiBr)**Koef_n[18]*(T_c/(T-T_0))**Koef_t[18]+\
        Koef_a[19]*X_mol_LiBr**Koef_m[19]*(0.4-X_mol_LiBr)**Koef_n[19]*(T_c/(T-T_0))**Koef_t[19]+\
        Koef_a[20]*X_mol_LiBr**Koef_m[20]*(0.4-X_mol_LiBr)**Koef_n[20]*(T_c/(T-T_0))**Koef_t[20]+\
        Koef_a[21]*X_mol_LiBr**Koef_m[21]*(0.4-X_mol_LiBr)**Koef_n[21]*(T_c/(T-T_0))**Koef_t[21]+\
        Koef_a[22]*X_mol_LiBr**Koef_m[22]*(0.4-X_mol_LiBr)**Koef_n[22]*(T_c/(T-T_0))**Koef_t[22]+\
        Koef_a[23]*X_mol_LiBr**Koef_m[23]*(0.4-X_mol_LiBr)**Koef_n[23]*(T_c/(T-T_0))**Koef_t[23]+\
        Koef_a[24]*X_mol_LiBr**Koef_m[24]*(0.4-X_mol_LiBr)**Koef_n[24]*(T_c/(T-T_0))**Koef_t[24]+\
        Koef_a[25]*X_mol_LiBr**Koef_m[25]*(0.4-X_mol_LiBr)**Koef_n[25]*(T_c/(T-T_0))**Koef_t[25]+\
        Koef_a[26]*X_mol_LiBr**Koef_m[26]*(0.4-X_mol_LiBr)**Koef_n[26]*(T_c/(T-T_0))**Koef_t[26]+\
        Koef_a[27]*X_mol_LiBr**Koef_m[27]*(0.4-X_mol_LiBr)**Koef_n[27]*(T_c/(T-T_0))**Koef_t[27]+\
        Koef_a[28]*X_mol_LiBr**Koef_m[28]*(0.4-X_mol_LiBr)**Koef_n[28]*(T_c/(T-T_0))**Koef_t[28])-s
    T = newton(fun,x0=300)
    return  T     


def  Calc_X_from_T_p_satLiBrSol_Patek(T,p):
    """
    Calc_X_from_T_p_LiBrSol_Patek
    Uses coefficients and formula obtained from Patek 2006
    Input:
        -   Temperature of LiBR solution                                [K]
        -   Pressure of saturated solution                              [Pa]
    Output:
        -   Molar concentration of LiBr in saturated solution           [-]
    """
    ## Constants
    T_c = 647.096              #[K]
    p_c = 22.064 * 10**6        #[Pa]
    # Table 4
    Koef_a = [  -2.41303    *   10**2,
                1.91750     *   10**7,
                -1.75521    *   10**8,
                3.25430     *   10**7,
                3.92571     *   10**2,
                -2.12626    *   10**3,
                1.85127     *   10**8,
                1.91216     *   10**3]
    Koef_t = [0,0,0,0,1,1,1,1]
    Koef_n = [0,5,6,3,0,2,6,0]
    Koef_m = [3,4,4,8,1,1,4,6]
    # Table 11
    Koef_beta = [1.0,1.5,3.0,3.5,4.0,7.5]
    Koef_alpha = [  -7.85951783,
                    1.84408259,
                    -11.7866497,
                    22.6807411,
                    -15.9618719,
                    1.80122502]
    ## Calculation
    # Calculate T_sat
    funT = lambda T:p_c*math.exp((T_c/T)*\
        (Koef_alpha[0]*(1-(T/T_c))**Koef_beta[0]+\
        Koef_alpha[1]*(1-(T/T_c))**Koef_beta[1]+\
        Koef_alpha[2]*(1-(T/T_c))**Koef_beta[2]+\
        Koef_alpha[3]*(1-(T/T_c))**Koef_beta[3]+\
        Koef_alpha[4]*(1-(T/T_c))**Koef_beta[4]+\
        Koef_alpha[5]*(1-(T/T_c))**Koef_beta[5]))-p
    T_sat = root_scalar(funT,bracket=[100,500]).root
    #T_sat = newton(funT,x0=300)
    T_R = T/T_c
    funX = lambda X: T-(Koef_a[0]*X**Koef_m[0]*(0.4-X)**Koef_n[0]*(T_R)**Koef_t[0]\
        +Koef_a[1]*X**Koef_m[1]*(0.4-X)**Koef_n[1]*(T_R)**Koef_t[1]\
        +Koef_a[2]*X**Koef_m[2]*(0.4-X)**Koef_n[2]*(T_R)**Koef_t[2]\
        +Koef_a[3]*X**Koef_m[3]*(0.4-X)**Koef_n[3]*(T_R)**Koef_t[3]\
        +Koef_a[4]*X**Koef_m[4]*(0.4-X)**Koef_n[4]*(T_R)**Koef_t[4]\
        +Koef_a[5]*X**Koef_m[5]*(0.4-X)**Koef_n[5]*(T_R)**Koef_t[5]\
        +Koef_a[6]*X**Koef_m[6]*(0.4-X)**Koef_n[6]*(T_R)**Koef_t[6]\
        +Koef_a[7]*X**Koef_m[7]*(0.4-X)**Koef_n[7]*(T_R)**Koef_t[7])-T_sat
    X = newton(funX,x0=0.3)
    #X = root_scalar(funX,bracket=[0.01,0.9]).root
    return  X 


def  Calc_state_SHEX_exit(T_sat, p_cond, h_des_in, m_rich, w_H2O_rich, stepsize):
    """
    Calc_state_SHEX_exit
    Input:
        -   Saturation temperature of Solution T                    [K]
        -   Condensation pressure                                   [Pa]
        -   Enthalpy of solution decorber inlet                     [J/kg]
        -   Rich solution mass flow rate                            [kg/s]
        -   Rich solution H2O concentration                         [-]
        -   Stepzise                                                [-]
    Output:
        -   Solution temperature desorber inlet                     [K]
    ---------------------------------------------------------------------- 
    Since the solution is heated above saturation, some refrigerant must be
    evaporated
    The solution temperature is increased until mass and energy equilibrium
    are reached
    """
    ## Definition of constants
    # ------------------------Necessary Constants---------------------------- #
    M_LiBr = 0.08685               #[kg/mol]
    M_H2O = 0.018015268            #[kg/mol]
    ## Calculation
    # saturation as starting point
    x_LiBr = Calc_X_from_T_p_satLiBrSol_Patek(T_sat,p_cond)
    h_mixture = Calc_h_from_T_X_LiBrSol_Patek(T_sat,x_LiBr)
    T_des_in = T_sat + stepsize
    # loop
    while  h_mixture<h_des_in: # Increase the temerature until energy equilibrium is achieved
        x_LiBr = Calc_X_from_T_p_satLiBrSol_Patek(T_des_in,p_cond)
        w_LiBr = x_LiBr*M_LiBr/ (x_LiBr*M_LiBr + (1-x_LiBr)*M_H2O)
        w_H2O = 1 - w_LiBr
        m_H2O = (w_H2O_rich - w_H2O)*m_rich
        m_sol = m_rich - m_H2O
        h_sol_mol = Calc_h_from_T_X_LiBrSol_Patek(T_des_in,x_LiBr)
        h_sol = h_sol_mol/(x_LiBr*M_LiBr + (1-x_LiBr)*M_H2O)
        h_H2O = CP.PropsSI('H','T',T_des_in,'P',p_cond,'Water')
        h_mixture = (m_sol*h_sol + m_H2O*h_H2O)/(m_sol+m_H2O)
        if  h_mixture<h_des_in:
            T_des_in = T_des_in + stepsize
    return  T_des_in 


def  Calc_state_valve_exit(T_sat, p_evap, h_abs_in, m_poor, w_H2O_poor, stepsize):
    """
    Calc_state_valve_exit
    Input:
        -   Saturation temperature of Solution T                    [K]
        -   Evaporation pressure                                    [Pa]
        -   Enthalpy of solution absorber inlet                     [J/kg]
        -   Poor solution mass flow rate                            [kg/s]
        -   Poor solution H2O concentration                         [-]
        -   Stepzise                                                [-]
    Output:
        -   Solution temperature absorber inlet                     [K]
    ---------------------------------------------------------------------- 
    Since the solution is depressurized beyond saturation some refrigerant
    has to be evaporated
    The solution temperature is decreases until mass and energy equilibrium
    are reached
    """
    ## Definition of constants
    # ------------------------Necessary Constants---------------------------- #
    M_LiBr = 0.08685               #[kg/mol]
    M_H2O = 0.018015268            #[kg/mol]
    ## Calculation
    # saturation as starting point
    x_LiBr = Calc_X_from_T_p_satLiBrSol_Patek(T_sat,p_evap)
    h_mixture = Calc_h_from_T_X_LiBrSol_Patek(T_sat,x_LiBr)
    T_abs_in = T_sat - stepsize
    # loop
    while  h_mixture>h_abs_in: # Decrease the temerature until energy equilibrium is achieved
        x_LiBr = Calc_X_from_T_p_satLiBrSol_Patek(T_abs_in,p_evap)
        w_LiBr = x_LiBr*M_LiBr/ (x_LiBr*M_LiBr + (1-x_LiBr)*M_H2O)
        w_H2O = 1 - w_LiBr
        m_H2O = (w_H2O_poor - w_H2O)*m_poor
        m_sol = m_poor - m_H2O
        h_sol_mol = Calc_h_from_T_X_LiBrSol_Patek(T_abs_in,x_LiBr)
        h_sol = h_sol_mol/(x_LiBr*M_LiBr + (1-x_LiBr)*M_H2O)
        h_H2O = CP.PropsSI('H','T',T_abs_in,'P',p_evap,'Water')
        h_mixture = (m_sol*h_sol + m_H2O*h_H2O)/(m_sol+m_H2O)
        if  h_mixture < h_abs_in:
            T_abs_in = T_abs_in - stepsize
    return  T_abs_in    


def crystallization_H2OLiBr(decider,value):
    """
    crystallization
    Calculates critical temperatures or concentrations for given value
    Use "T" if Input is temperature as decider value
    Use "w" if Input is mass fraction of LiBr in solution as decider value
    """
    # # Computation
    T_cr = 0
    w_cr = 0
    if decider == "w":
        # Parameters (from Albers report 2019, Boryta 1970)
        parameter_T = [ 42.90198341384762,
                        34.67510890651030,
                        31.30778644395644,
                        2.99859601946791,
                        -19.36781324384540,
                        -4.88856108511827,
                        4.61433775768846,
                        1.80636830673333   ]
        # Define working variable
        w = value
        w_r = (w-0.64794)/0.044858;
        # Calculate critical values
        for i in range(8):
            T_cr = T_cr + parameter_T[i]*(w_r**(i))
        CritValue = T_cr
    elif decider == "T":
        # Parameters (from Albers report 2019, Boryta 1970)
        parameter_w = [ 0.66136507494441,
                        0.02262634534253,
                        -0.02216522722755,
                        0.05134156572205,
                        0.00034455919818,
                        -0.03628931060739,
                        0.00252166562759,
                        0.00796985214167    ]
        # Define working variable
        T = value - 273.15
        T_r = (T-54.793)/33.111
        # Calculate critical values
        for i in range(8):
            w_cr = w_cr + parameter_w[i]*(T_r**(i))
        CritValue = w_cr
    else:
        print("Crystallization function failed - Decider value not defined - use T or w")
        return
    return CritValue


def checkForViolation_H2OLiBr(w,T,position):
    """
    checkForViolation
    Checks for violation of accepted values for solution state functions
    Use "T" if Input is temperature as decider value
    Use "w" if Input is mass fraction of LiBr in solution as decider value
    """
    # # Computation
    w_cr = crystallization_H2OLiBr("T",T)
    T_cr = crystallization_H2OLiBr("w",w)
    # Check for Patek boundaries (273K - 500K)
    if T>500:
        print(f"Temperature of LiBr solution ({T} K) above Patek max. tempreture of 500 K at point: {position}")
        return
    elif T<=273.15:
        print(f"Temperature of LiBr solution ( {T} K) below Patek min. tempreture of 273.15 K at point: {position}")
        return 
        
    
    # Check for Boryta boundaries (273K - 374K)
    if T>374:
        print(f"Crystallization can not be checked because Temperature (",T,"K) of LiBr solution above Boryta max. tempreture of 374 K at point: ",position)
    elif T<=273.15:
        print(f"Crystallization can not be checked because Temperature of LiBr solution (",T,"K) below min. tempreture of 273.15 K at point: ",position)
    
    if T<=374 and T>=273.15 and w>0.57:
        if w>w_cr:
            print(f"Crystallization Error! Mass fraction of LiBr in solution({w}) above max. value of {w_cr} at position: {position}")
            return
        elif T<T_cr:
            print(f"Crystallization Error! Temperature of LiBr in solution ({T}K) below min. value of {T_cr}K at position: {position}")
            return
        
    return

def find_starting_value(df, variable, prop1, prop2, value1, value2):
    """ 
        df = dataframe containing the data for comparison;
        variable = str of property for which the start value is searched for
        prop1, prop2 = str of the properties provided to find the starting value
        value1, value2 = float of the values of prop1 and prop2

        Is used to finde the closest starting value for solving the inverse functions
        using the newton method
    """
    values = df[prop1].unique()
    Delta_v1, v1_close, Delta_v2, v2_close = 10000, 0, 10000, 0
    for v1 in values:
        if Delta_v1 > abs(value1-v1):
            Delta_v1 =  abs(value1-v1)
            v1_close = v1
    
    df = df[df[prop1]==v1_close]
    values = df[prop2].unique()
    for v2 in values:
        if Delta_v2 > abs(value2-v2):
            Delta_v2 = abs(value2-v2)
            v2_close = v2

    df = df[df[prop2] == v2_close]

    return df.iloc[0][variable]


def Calc_X_from_T_cp_LiBrSol(T,cp, X_start=0):
    """ 
        Calculation of:
            X       Molar Concentration X of LiBr in Solution                   [-]                [-]

        Based on:
            T       Temperature of the solution                           [K]
            cp      Molar heat capacity                                   [mol/m**3]

        Solves the inverse function numerical, by using the newton method     
    """
    
    
    X_start = find_starting_value(DFcal_LiBrSol, "X","T [K]", "cp [J/molK]", T, cp)
    problem = lambda X, T, cp: (Calc_cp_from_T_X_LiBrSol_Patek(T,X)-cp)**2

    X_fit = newton(problem, x0 = X_start, args=(T, cp))

    return X_fit

def Calc_T_from_X_cp_LiBrSol(X,cp, T_start=10_000):
    """ 
        Calculation of:
            T       Temperature of the solution                           [K]               

        Based on:
            X       Molar Concentration X of LiBr in Solution                   [-]                [-]
            cp      Molar heat capacity                                   [mol/m**3]

        Solves the inverse function numerical, by using the newton method     
    """
    T_start = find_starting_value(DFcal_LiBrSol, "T [K]" ,"X", "cp [J/molK]", X, cp)
    problem = lambda T, X, cp: (Calc_cp_from_T_X_LiBrSol_Patek(T,X)-cp)**2

    T_fit = newton(problem, x0 = T_start, args=(X, cp))

    return T_fit

def Calc_X_from_T_h_LiBrSol(T,h, X_start=10_000):
    """ 
        Calculation of:
            X       Molar Concentration X of LiBr in Solution                   [-]                [-]              

        Based on:
            T       Temperature of the solution                           [K]
            h       Molar enthalpy h of solution                          [J/mol]   

        Solves the inverse function numerical, by using the newton method     
    """
    X_start = find_starting_value(DFcal_LiBrSol, "X","T [K]", "h [J/mol]", T, h)
    problem = lambda X, T, h: (Calc_h_from_T_X_LiBrSol_Patek(T,X)-h)**2

    X_fit = newton(problem, x0 = X_start, args=(T, h))

    return X_fit

def Calc_T_from_X_h_LiBrSol(X,h, T_start=10_000):
    """ 
        Calculation of:
            T       Temperature of the solution                           [K]         

        Based on:
            X       Molar Concentration X of LiBr in Solution                   [-]                [-]              
            h       Molar enthalpy h of solution                          [J/mol]   

        Solves the inverse function numerical, by using the newton method     
    """
    T_start = find_starting_value(DFcal_LiBrSol,  "T [K]" ,"X", "h [J/mol]", X, h)
    problem = lambda T, X, h: (Calc_h_from_T_X_LiBrSol_Patek(T,X)-h)**2

    T_fit = newton(problem, x0 = T_start, args=(X, h))

    return T_fit

def Calc_X_from_T_p_LiBrSol(T, p, X_start=10_000):
    """ 
        Calculation of:
            X       Molar Concentration X of LiBr in Solution                   [-]                [-]       

        Based on:
            T       Temperature of the solution                           [K]               
            p       Pressure of saturated solution                        [Pa]

        Solves the inverse function numerical, by using the newton method     
    """
    
    X_start = find_starting_value(DFcal_LiBrSol, "X","T [K]", "p [Pa]", T, p)
    problem = lambda X, T, p: (Calc_p_from_T_X_LiBrSol_Patek(T,X)-p)**2

    X_fit = newton(problem, x0 = X_start, args=(T, p))

    return X_fit

def Calc_T_from_X_p_LiBrSol(X,p, T_start=10_000):
    """ 
        Calculation of:
            
            T       Temperature of the solution                           [K]
        Based on:
            X       Molar Concentration X of LiBr in Solution                   [-]                [-]                      
            p       Pressure of saturated solution                        [Pa]

        Solves the inverse function numerical, by using the newton method     
    """
    T_start = find_starting_value(DFcal_LiBrSol,  "T [K]" ,"X", "p [Pa]", X, p)
    problem = lambda T, X, p: (Calc_p_from_T_X_LiBrSol_Patek(T,X)-p)**2

    T_fit = newton(problem, x0 = T_start, args=(X, p))

    return T_fit

def Calc_X_from_T_rho_LiBrSol(T,rho, X_start=10_000):
    """ 
    Calculation of:
        X       Molar Concentration X of LiBr in Solution             [-]

    Based on:
        T       Temperature of the solution                           [K]                      
        rho     Molar density rho                                     [mol/m**3]

    Solves the inverse function numerical, by using the newton method    
    """

    X_start = find_starting_value(DFcal_LiBrSol,"X", "T [K]", "rho [mol/m3]", T, rho)
    problem = lambda X, T, rho: (Calc_rho_from_T_X_LiBrSol_Patek(T,X)-rho)**2

    X_fit = newton(problem, x0 = X_start, args=(T, rho))

    return X_fit

def Calc_T_from_X_rho_LiBrSol(X,rho, T_start=10_000):
    """ 
    Calculation of:           
        T       Temperature of the solution                           [K]     

    Based on:
        X       Molar Concentration X of LiBr in Solution             [-]                      
        rho     Molar density rho                                     [mol/m**3]

    Solves the inverse function numerical, by using the newton method    
    """
    T_start = find_starting_value(DFcal_LiBrSol,  "T [K]" ,"X", "rho [mol/m3]", X, rho)
    problem = lambda T, X, rho: (Calc_rho_from_T_X_LiBrSol_Patek(T,X)-rho)**2

    T_fit = newton(problem, x0 = T_start, args=(X, rho))

    return T_fit

def Calc_X_from_T_s_LiBrSol(T,s, X_start=10_000):
    """ 
    Calculation of:              
        X       Molar Concentration X of LiBr in Solution             [-]     

    Based on:
        T       Temperature of the solution                           [K]                    
        s       Molar enthalpy of solution                            [J/molK]

    Solves the inverse function numerical, by using the newton method    
    """
    X_start = find_starting_value(DFcal_LiBrSol,"X", "T [K]", "s [J/molK]", T, s)
    problem = lambda X, T, s: (Calc_s_from_T_X_LiBrSol_Patek(T,X)-s)**2

    X_fit = newton(problem, x0 = X_start, args=(T, s))

    return X_fit

def Calc_T_from_X_s_LiBrSol(X,s, T_start=10_000):
    """ 
    Calculation of:
        T       Temperature of the solution                           [K]             
        
    Based on:
        X       Molar Concentration X of LiBr in Solution             [-]                         
        s       Molar enthalpy of solution                            [J/molK]

    Solves the inverse function numerical, by using the newton method    
    """
    T_start = find_starting_value(DFcal_LiBrSol,  "T [K]" ,"X", "s [J/molK]", X, s)
    problem = lambda T, X, s: (Calc_s_from_T_X_LiBrSol_Patek(T,X)-s)**2

    T_fit = newton(problem, x0 = T_start, args=(X, s))

    return T_fit
