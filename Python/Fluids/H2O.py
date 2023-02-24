# class H2O():
#     """
#     Class containing the calculations for H2O
#     """

#     def __init__( math, newton, pandas):
#         """
#         Initialisation of the used moduels. 
#         The needen models are:
#             - math
#             - the funciton newton o the scipy package
#             - pandas
#         To improve performance, the moduels are imported once and will be provided to the class at initialisation.
        
#         In Addition the look up table of starting values for solving the invers functions is imported 
#         and stored in a data frame.
#         """
#         math = math
#         newton = newton
#         pd = pandas
#         DFcal = pd.read_csv("01_Class_approach\H2O\H2O_cal.csv", sep = ";")
import pandas as pd
from scipy.optimize import newton, fsolve
import math
import numpy
import CoolProp.CoolProp as CP
global DFcal_H2O 
DFcal_H2O = pd.read_csv("Python\Fluids\H2O_cal.csv", sep = ";")


def Calc_cp_from_T_p_H2O_Patek(T,p):
    """T [K], p [Pa]
    # ---------------------------------------------------------------------- #
    # Calc_cp_from_T_p_H2O_Patek
    # Uses coefficients and formula obtained from Patek 2009
    # ---------------------------------------------------------------------- #
    # {
    # Author  : Ludwig Irrgang
    # Date    : 11.10.2022
    # Copyright information:
    # Ludwig Irrgang
    # Lehrstuhl für Energiesysteme
    # TUM School of Engineering and Design
    # Technische Universität München
    # Boltzmannstr. 15 
    # 85748 Garching b. München
    # ludwig.irrgang@tum.de
    # }
    # ---------------------------------------------------------------------- #
    # Input:
    #       -   Termperature of Solution T                                  [K]
    #       -   pressure of the steam p                                     [Pa]
    # Output:
    #       -   molar isobaric heat capacity                                [J/molK]
    # ---------------------------------------------------------------------- #
    """
    ## Constants
    T_c = 647.096           # [K]
    p_c = 22.064 * 10**6    # [Pa]
    R = 8.314371            # [J/molK]
    M_H20 = 0.018015268     # [kg/mol]
    
    ## Table 1
    Koeff_a= [  4.03740,
                6.66860,
                -8.78474,
                -2.67859*10**-2,
                -8.24016*10**-1,
                7.39848*10**-1,
                -2.73558*10**-1,
                2.51325*10**-1,
                -8.64644*10**-2,
                4.04713*10**-2,
                -2.07714,
                2.48941,
                -1.26382,
                6.46755,
                -4.63404,
                -3.59518,
                1.71653]
    Koeff_m = [-1, 0,0,0,1,1,1,2,2,2,3,3,3,4,4,5,5]
    Koeff_n = [0,1,0,-4,4,5,6,2,12,13,12,15,16,15,16,15,16]

    ## Calculation
    # Calculation of tau an pi
    tau = T_c/T
    pi = p/p_c

    sum = 0
    for i in range(1,17,1):
        sum = sum + Koeff_n[i]*(Koeff_n[i]-1)*Koeff_a[i]*tau**Koeff_n[i]*pi**Koeff_m[i]
    
    c_p = R*(Koeff_a[0]-sum)

    return c_p


def Calc_g_from_T_p_H2O_Patek(T,p):
    """T [K], p [Pa]
    # ---------------------------------------------------------------------- #
    # Calc_g_from_T_p_H2O_Patek
    # Uses coefficients and formula obtained from Patek 2009
    # ---------------------------------------------------------------------- #
    # {
    # Author  : Ludwig Irrgang
    # Date    : 11.10.2022
    # Copyright information:
    # Ludwig Irrgang
    # Lehrstuhl für Energiesysteme
    # TUM School of Engineering and Design
    # Technische Universität München
    # Boltzmannstr. 15 
    # 85748 Garching b. München
    # ludwig.irrgang@tum.de
    # }
    # ---------------------------------------------------------------------- #
    # Input:
    #       -   Termperature of Solution T                                  [K]
    #       -   pressure of the steam p                                     [Pa]
    # Output:
    #       -   molar Gibbs free energy                                     [J/mol]
    # ---------------------------------------------------------------------- #
    """
    ## Constants
    T_c = 647.096           # [K]
    p_c = 22.064 * 10**6    # [Pa]
    R = 8.314371            # [J/molK]
    M_H20 = 0.018015268     # [kg/mol]
    
    ## Table 1
    Koeff_a= [  4.03740,
                6.66860,
                -8.78474,
                -2.67859*10**-2,
                -8.24016*10**-1,
                7.39848*10**-1,
                -2.73558*10**-1,
                2.51325*10**-1,
                -8.64644*10**-2,
                4.04713*10**-2,
                -2.07714,
                2.48941,
                -1.26382,
                6.46755,
                -4.63404,
                -3.59518,
                1.71653]
    Koeff_m = [-1, 0,0,0,1,1,1,2,2,2,3,3,3,4,4,5,5]
    Koeff_n = [0,1,0,-4,4,5,6,2,12,13,12,15,16,15,16,15,16]

    ## Calculation
    # Calculation of tau an pi
    tau = T_c/T
    pi = p/p_c

    sum = 0
    for i in range(1,17,1):
        sum = sum + Koeff_a[i]*tau**Koeff_n[i]*pi**Koeff_m[i]
    
    g = R*T*(math.log(pi)+Koeff_a[0]*math.log(tau)+sum)

    return g

def Calc_h_from_T_p_H2O_Patek(T,p):
    """
    T [K], p [Pa]
    # ---------------------------------------------------------------------- #
    # Calc_h_from_T_p_H2O_Patek
    # Uses coefficients and formula obtained from Patek 2009
    # ---------------------------------------------------------------------- #
    # {
    # Author  : Ludwig Irrgang
    # Date    : 11.10.2022
    # Copyright information:
    # Ludwig Irrgang
    # Lehrstuhl für Energiesysteme
    # TUM School of Engineering and Design
    # Technische Universität München
    # Boltzmannstr. 15 
    # 85748 Garching b. München
    # ludwig.irrgang@tum.de
    # }
    # ---------------------------------------------------------------------- #
    # Input:
    #       -   Termperature of Solution T                                  [K]
    #       -   pressure of the steam p                                     [Pa]
    # Output:
    #       -   molar enthalpy                                              [J/mol]
    # ---------------------------------------------------------------------- #
    """

    ## Constants
    T_c = 647.096           # [K]
    p_c = 22.064 * 10**6    # [Pa]
    R = 8.314371            # [J/molK]
    M_H20 = 0.018015268     # [kg/mol]
    
    ## Table 1
    Koeff_a= [  4.03740,
                6.66860,
                -8.78474,
                -2.67859*10**-2,
                -8.24016*10**-1,
                7.39848*10**-1,
                -2.73558*10**-1,
                2.51325*10**-1,
                -8.64644*10**-2,
                4.04713*10**-2,
                -2.07714,
                2.48941,
                -1.26382,
                6.46755,
                -4.63404,
                -3.59518,
                1.71653]
    Koeff_m = [-1, 0,0,0,1,1,1,2,2,2,3,3,3,4,4,5,5]
    Koeff_n = [0,1,0,-4,4,5,6,2,12,13,12,15,16,15,16,15,16]

    ## Calculation
    # Calculation of tau an pi
    tau = T_c/T
    pi = p/p_c

    sum = 0
    for i in range(1,17,1):
        sum = sum + Koeff_n[i]*Koeff_a[i]*tau**Koeff_n[i]*pi**Koeff_m[i]
    
    h = R*T*(Koeff_a[0]+sum)

    return h

def Calc_s_from_T_p_H2O_Patek(T,p):
    """
    T [K], p [Pa]
    # ---------------------------------------------------------------------- #
    # Calc_s_from_T_p_H2O_Patek
    # Uses coefficients and formula obtained from Patek 2009
    # ---------------------------------------------------------------------- #
    # {
    # Author  : Ludwig Irrgang
    # Date    : 11.10.2022
    # Copyright information:
    # Ludwig Irrgang
    # Lehrstuhl für Energiesysteme
    # TUM School of Engineering and Design
    # Technische Universität München
    # Boltzmannstr. 15 
    # 85748 Garching b. München
    # ludwig.irrgang@tum.de
    # }
    # ---------------------------------------------------------------------- #
    # Input:
    #       -   Termperature of Solution T                                  [K]
    #       -   pressure of the steam p                                     [Pa]
    # Output:
    #       -   molar entropy                                               [J/molK]
    # ---------------------------------------------------------------------- #
    """
    ## Constants
    T_c = 647.096           # [K]
    p_c = 22.064 * 10**6    # [Pa]
    R = 8.314371            # [J/molK]
    M_H20 = 0.018015268     # [kg/mol]
    
    ## Table 1
    Koeff_a= [  4.03740,
                6.66860,
                -8.78474,
                -2.67859*10**-2,
                -8.24016*10**-1,
                7.39848*10**-1,
                -2.73558*10**-1,
                2.51325*10**-1,
                -8.64644*10**-2,
                4.04713*10**-2,
                -2.07714,
                2.48941,
                -1.26382,
                6.46755,
                -4.63404,
                -3.59518,
                1.71653]
    Koeff_m = [-1, 0,0,0,1,1,1,2,2,2,3,3,3,4,4,5,5]
    Koeff_n = [0,1,0,-4,4,5,6,2,12,13,12,15,16,15,16,15,16]

    ## Calculation
    # Calculation of tau an pi
    tau = T_c/T
    pi = p/p_c

    sum = 0
    for i in range(1,17,1):
        sum = sum + (Koeff_n[i]-1)*Koeff_a[i]*tau**Koeff_n[i]*pi**Koeff_m[i]
    
    s = R*(-math.log(pi)+Koeff_a[0]*(1-math.log(tau))+sum)

    return s

def Calc_v_from_T_p_H2O_Patek(T,p):
    """T [K], p [Pa]
    # ---------------------------------------------------------------------- #
    # Calc_v_from_T_p_H2O_Patek
    # Uses coefficients and formula obtained from Patek 2009
    # ---------------------------------------------------------------------- #
    # {
    # Author  : Ludwig Irrgang
    # Date    : 11.10.2022
    # Copyright information:
    # Ludwig Irrgang
    # Lehrstuhl für Energiesysteme
    # TUM School of Engineering and Design
    # Technische Universität München
    # Boltzmannstr. 15 
    # 85748 Garching b. München
    # ludwig.irrgang@tum.de
    # }
    # ---------------------------------------------------------------------- #
    # Input:
    #       -   Termperature of Solution T                                  [K]
    #       -   pressure of the steam p                                     [Pa]
    # Output:
    #       -   molar volume                                                [m3/mol]
    # ---------------------------------------------------------------------- #
    """
    ## Constants
    T_c = 647.096           # [K]
    p_c = 22.064 * 10**6    # [Pa]
    R = 8.314371            # [J/molK]
    M_H20 = 0.018015268     # [kg/mol]
    
    ## Table 1
    Koeff_a= [  4.03740,
                6.66860,
                -8.78474,
                -2.67859*10**-2,
                -8.24016*10**-1,
                7.39848*10**-1,
                -2.73558*10**-1,
                2.51325*10**-1,
                -8.64644*10**-2,
                4.04713*10**-2,
                -2.07714,
                2.48941,
                -1.26382,
                6.46755,
                -4.63404,
                -3.59518,
                1.71653]
    Koeff_m = [-1, 0,0,0,1,1,1,2,2,2,3,3,3,4,4,5,5]
    Koeff_n = [0,1,0,-4,4,5,6,2,12,13,12,15,16,15,16,15,16]

    ## Calculation
    # Calculation of tau an pi
    tau = T_c/T
    pi = p/p_c

    sum = 0
    for i in range(1,17,1):
        sum = sum + Koeff_m[i]*Koeff_a[i]*tau**Koeff_n[i]*pi**Koeff_m[i]
    
    v = (R*T/p)*(1+sum)

    return v

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

def Calc_p_from_T_cp_H2O(T,cp, p_start=10_000):
    """ 
    Calculation of:
        p       pressure of the steam                           [Pa]             
        
    Based on:
        T       Temperature of solution                         [K]                         
        cp      molar isobaric heat capacity                    [J/molK]

    Solves the inverse function numerical, by using the newton method    
    """
    p_start = find_starting_value(DFcal_H2O, "p [Pa]", "T [K]", "cp [J/molK]", T, cp)

    problem = lambda p, T, cp: (Calc_cp_from_T_p_H2O_Patek(T,p)-cp)**2

    p_fit = newton(problem, x0 = p_start, args=(T, cp))

    return p_fit

def Calc_T_from_p_cp_H2O(p,cp, T_start=10_000):
    """ 
    Calculation of:                        
        T       Temperature of solution                         [K]

    Based on:
        p       pressure of the steam                           [Pa]                          
        cp      molar isobaric heat capacity                    [J/molK]

    Solves the inverse function numerical, by using the newton method    
    """
    T_start = find_starting_value(DFcal_H2O, "T [K]", "p [Pa]", "cp [J/molK]", p, cp)
    problem = lambda T, p, cp: (Calc_cp_from_T_p_H2O_Patek(T,p)-cp)**2

    T_fit = newton(problem, x0 = T_start, args=(p, cp))

    return T_fit

def Calc_p_from_T_g_H2O(T,g, p_start=10_000):
    """ 
    Calculation of:
        p       pressure of the steam                           [Pa]             
        
    Based on:
        T       Temperature of solution                         [K]                         
        g       molar Gibbs free energy                         [J/mol]

    Solves the inverse function numerical, by using the newton method    
    """
    p_start = find_starting_value(DFcal_H2O, "p [Pa]", "T [K]", "g [J/mol]", T,g)

    problem = lambda p, T, g: (Calc_g_from_T_p_H2O_Patek(T,p)-g)**2

    p_fit = newton(problem, x0 = p_start, args=(T, g))

    return p_fit

def Calc_T_from_p_g_H2O(p,g, T_start=10_000):
    """ 
    Calculation of:                       
        T       Temperature of solution                         [K]

    Based on:
        p       pressure of the steam                           [Pa]                         
        g       molar Gibbs free energy                         [J/mol]

    Solves the inverse function numerical, by using the newton method    
    """
    T_start = find_starting_value(DFcal_H2O, "T [K]", "p [Pa]", "g [J/mol]", p,g)
    problem = lambda T, p, g: (Calc_g_from_T_p_H2O_Patek(T,p)-g)**2

    T_fit = newton(problem, x0 = T_start, args=(p, g))

    return T_fit

def Calc_p_from_T_h_H2O(T,h, p_start=10_000):
    """ 
    Calculation of:
        p       pressure of the steam                           [Pa]             
        
    Based on:
        T       Temperature of solution                         [K]                         
        h       molar enthalpy                                  [J/mol]

    Solves the inverse function numerical, by using the newton method    
    """
    p_start = find_starting_value(DFcal_H2O, "p [Pa]", "T [K]", "h [J/mol]", T, h)
    problem = lambda p, T, h: (Calc_h_from_T_p_H2O_Patek(T,p)-h)**2

    p_fit = newton(problem, x0 = p_start, args=(T, h))

    return p_fit

def Calc_T_from_p_h_H2O(p,h, T_start=10_000):
    """ 
    Calculation of:                         
        T       Temperature of solution                         [K]

    Based on:
        p       pressure of the steam                           [Pa]                        
        h       molar enthalpy                                  [J/mol]

    Solves the inverse function numerical, by using the newton method    
    """
    T_start = find_starting_value(DFcal_H2O, "T [K]", "p [Pa]", "h [J/mol]", p , h)
    problem = lambda T, p, h: (Calc_h_from_T_p_H2O_Patek(T,p)-h)**2

    T_fit = newton(problem, x0 = T_start, args=(p, h))

    return T_fit

def Calc_p_from_T_s_H2O(T,s, p_start=10_000):
    """ 
    Calculation of:
        p       pressure of the steam                           [Pa]             
        
    Based on:
        T       Temperature of solution                         [K]                         
        s       molar entropy                                   [J/molK]

    Solves the inverse function numerical, by using the newton method    
    """
    p_start = find_starting_value(DFcal_H2O, "p [Pa]", "T [K]", "s [J/molK]", T, s)
    problem = lambda p, T, s: (Calc_s_from_T_p_H2O_Patek(T,p)-s)**2

    p_fit = newton(problem, x0 = p_start, args=(T, s))

    return p_fit

def Calc_T_from_p_s_H2O(p,s, T_start=10_000):
    """ 
    Calculation of:                       
        T       Temperature of solution                         [K]

    Based on:
        p       pressure of the steam                           [Pa]                           
        s       molar entropy                                   [J/molK]

    Solves the inverse function numerical, by using the newton method    
    """
    T_start = find_starting_value(DFcal_H2O, "T [K]", "p [Pa]", "s [J/molK]", p, s)
    problem = lambda T, p, s: (Calc_s_from_T_p_H2O_Patek(T,p)-s)**2

    T_fit = newton(problem, x0 = T_start, args=(p, s))

    return T_fit

def Calc_p_from_T_v_H2O(T,v, p_start=10_000):
    """ 
    Calculation of:
        p       pressure of the steam                           [Pa]             
        
    Based on:
        T       Temperature of solution                         [K]                         
        v       molar volume                                    [m3/mol]

    Solves the inverse function numerical, by using the newton method    
    """
    p_start = find_starting_value(DFcal_H2O, "p [Pa]", "T [K]", "v [m3/mol]", T, v)
    problem = lambda p, T, v: (Calc_v_from_T_p_H2O_Patek(T,p)-v)**2

    p_fit = newton(problem, x0 = p_start, args=(T, v))

    return p_fit

def Calc_T_from_p_v_H2O(p,v, T_start=10_000):
    """ 
    Calculation of:                       
        T       Temperature of solution                         [K]

    Based on:
        p       pressure of the steam                           [Pa]                           
        v       molar volume                                    [m3/mol]
    Solves the inverse function numerical, by using the newton method    
    """
    T_start = find_starting_value(DFcal_H2O, "T [K]", "p [Pa]", "v [m3/mol]", p, v)
    problem = lambda T, p, v: (Calc_cp_from_T_p_H2O_Patek(T,p)-v)**2

    T_fit = newton(problem, x0 = T_start, args=(p, v))

    return T_fit

