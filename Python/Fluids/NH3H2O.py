import math
import sys
from scipy.optimize import newton, ridder, bisect
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import os

os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])

def NH3_conversion_w_X(direction, a):
    """
    input:
        direction: "X" or "Y" for conversion to mol fraction
                   "w" for conversion to mass fraction
    Converts mole fraction to mass fraction and vice versa
    """
    M_w = 18.015        # [g/mol]
    M_a = 17.031        # [g/mol]
    if a == 0:
        b = 0
    elif direction == "Y" or direction =="X":
        b = 1/((1/M_w)*(M_a/a - M_a)+1)
    elif direction =="w":
        b = (a*M_a)/(a*M_a + (1-a)*M_w)

    return b


def Calc_T_from_p_X_NH3H2OSol_Patek( p,X_mol_NH3):
    """
    % ---------------------------------------------------------------------- %
    % Calc_T_from_p_X_NH3H2OSol_Patek
    % Uses coefficients and formula obtained from Patek 1995
    % ---------------------------------------------------------------------- %
    %{
    Author  : Ludwig Irrgang
    Date    : 25.06.2022
    Copyright information:
    Ludwig Irrgang
    Lehrstuhl für Energiesysteme
    TUM School of Engineering and Design
    Technische Universität München
    Boltzmannstr. 15 
    85748 Garching b. München
    ludwig.irrgang@tum.de
    %}
    % ---------------------------------------------------------------------- %
    % Input:
    %       -   Pressure                                                   [Pa]
    %       -   Molar Concentration X of NH3 in Solution                   [-]
    % Output:
    %       -   Temperature of NH3 - H2O solution                          [K]
    % ---------------------------------------------------------------------- %
    """
    ## Unit Conversion
    p = p*10**(-6) # from Pa to MPa for calculation with formulas from Patek 1995
    ## Constants
    p_0 = 2                   # [MPa]
    T_0 = 100                   #[K]
    # Table 1
    Koef_a=[ 0.322302 * 10**1,
                -0.384206*10**0,
                0.460965 * 10 **-1,
                -0.378945 * 10**-2,
                0.135610 * 10**-3,
                0.487755 * 10**0,
                -0.120108 * 10 **0,
                0.106154 * 10**-1,
                -0.533589 * 10**-3,
                0.785041 * 10**1,
                -0.115941 *10**2,
                -0.523150*10**-1,
                0.489596 * 10**1,
                0.421059 *10**-1
    ]
    Koef_n = [0,1,2,3,4,0,1,2,3,0,0,1,0,1]
    Koef_m = [0,0,0,0,0,1,1,1,2,4,5,5,6,13]    

    ## Calculation
    sum = 0
    for i in range(14):
        sum = sum + (Koef_a[i]*(1-X_mol_NH3)**Koef_m[i])*( math.log(p_0/p))**Koef_n[i]

    
    T = T_0 * sum
    
    return T

def Calc_T_from_p_Y_NH3H2OSol_Patek( p,Y_mol_NH3):
    """
    % ---------------------------------------------------------------------- %
    % Calc_T_from_p_Y_NH3H2OSol_Patek
    % Uses coefficients and formula obtained from Patek 1995
    % ---------------------------------------------------------------------- %
    %{
    Author  : Ludwig Irrgang
    Date    : 25.06.2022
    Copyright information:
    Ludwig Irrgang
    Lehrstuhl für Energiesysteme
    TUM School of Engineering and Design
    Technische Universität München
    Boltzmannstr. 15 
    85748 Garching b. München
    ludwig.irrgang@tum.de
    %}
    % ---------------------------------------------------------------------- %
    % Input:
    %       -   Pressure                                                   [Pa]
    %       -   Molar Concentration Y of NH3 in the gas                    [-]
    % Output:
    %       -   Temperature of NH3 - H2O solution                          [K]
    % ---------------------------------------------------------------------- %
    """
    ## Unit Conversion
    p = p*10**(-6) # from Pa to MPa for calculation with formulas from Patek 1995
    ## Constants
    p_0 = 2                   # [MPa]
    T_0 = 100                   #[K]
    # Table 2
    Koef_a=[ 0.324004 * 10**1,
                -0.395920*10**0,
                0.435624 * 10 **-1,
                -0.218943 * 10**-2,
                -0.143526 * 10**1,
                0.105256 * 10**1,
                -0.719281 * 10 **-1,
                0.122362 * 10**2,
                -0.224368 * 10**1,
                -0.201780 * 10**2,
                0.110834 *10**1,
                0.145399*10**2,
                0.644312 * 10**0,
                -0.221246 *10**1,
                -0.756266*10**0,
                -0.135529 * 10**1,
                0.183541*10**0
    ]
    Koef_n = [0,1,2,3,0,1,2,0,1,0,1,0,2,0,2,0,2]
    Koef_m = [0,0,0,0,1,1,1,2,2,3,3,4,4,5,5,6,7]    

    ## Calculation
    sum = 0
    for i in range(17):
        sum = sum + (Koef_a[i]*(1-Y_mol_NH3)**(0.25*Koef_m[i]))*( math.log(p_0/p))**Koef_n[i]
    
    T = T_0 * sum
    
    return T

def Calc_Y_from_p_X_NH3H2OSol_Patek( p,X_mol_NH3):
    """
    % ---------------------------------------------------------------------- %
    % Calc_T_from_p_Y_NH3H2OSol_Patek
    % Uses coefficients and formula obtained from Patek 1995
    % ---------------------------------------------------------------------- %
    %{
    Author  : Ludwig Irrgang
    Date    : 25.06.2022
    Copyright information:
    Ludwig Irrgang
    Lehrstuhl für Energiesysteme
    TUM School of Engineering and Design
    Technische Universität München
    Boltzmannstr. 15 
    85748 Garching b. München
    ludwig.irrgang@tum.de
    %}
    % ---------------------------------------------------------------------- %
    % Input:
    %       -   Pressure                                                   [Pa]
    %       -   Molar Concentration X of NH3 in the solution               [-]
    % Output:
    %       -   Molar Concentration Y of NH3 in the solution               [-]
    % ---------------------------------------------------------------------- %
    """
    ## Unit Conversion
    p = p*10**(-6)  # from Pa to MPa for calculation with formulas from Patek 1995
    ## Constants
    p_0 = 2                   # [MPa]
    # Table 3
    Koef_a=[ 1.98022017*10,
                -1.18092669 * 10,
                2.77479980 * 10,
                -2.88634277 * 10,
                -5.91616608 * 10,
                5.78091305 *10**2,
                -6.21736743 *10**0,
                -3.42198402*10**3,
                1.19403127 * 10**4,
                -2.45413777 * 10**4,
                2.91591865*10**4,
                -1.84782290*10**4,
                2.34819434 *10,
                4.80310617*10**3
    ]
    Koef_n = [0,1,6,7,0,1,2,2,3,4,5,6,7,7]
    Koef_m = [0,0,0,0,1,2,2,3,4,5,6,7,7,8]    

    ## Calculation
    if X_mol_NH3==1:
        Y_mol_NH3 = 1
    else:
        sum = 0
        for i in range(14):
            sum = sum + (Koef_a[i]*(p/p_0)**Koef_m[i]*X_mol_NH3**(Koef_n[i]/3))

        Y_mol_NH3 = 1-  math.exp( math.log(1-X_mol_NH3)*sum)
    
    return Y_mol_NH3

def Calc_h_liquid_from_T_X_NH3H2O_Patek( T, X_mol_NH3):
    """
    % ---------------------------------------------------------------------- %
    % Calc_h_from_T_X_NH3H2O_Patek
    % Uses coefficients and formula obtained from Patek 1995
    % ---------------------------------------------------------------------- %
    %{
    Author  : Ludwig Irrgang
    Date    : 25.06.2022
    Copyright information:
    Ludwig Irrgang
    Lehrstuhl für Energiesysteme
    TUM School of Engineering and Design
    Technische Universität München
    Boltzmannstr. 15 
    85748 Garching b. München
    ludwig.irrgang@tum.de
    %}
    % ---------------------------------------------------------------------- %
    % Input:
    %       -   Termperature of Solution T                              [K]
    %       -   Molar Concentration X of NH3 in Solution                [-]
    % Output:
    %       -   specific enthalpy h of solution                         [J/kg]
    % ---------------------------------------------------------------------- %
    """
    ## Constants
    h_0 = 100                   # [kJ/kg]
    T_0 = 273.16                   #[K]
    # Table 4
    Koef_a=[  -0.761080 * 10**1,
                0.256905 * 10**2,
                -0.247092 * 10**3,
                0.325952 * 10**3,
                -0.158854 * 10**3,
                0.619084 * 10**2,
                0.114314 * 10**2,
                0.118157 * 10**1,
                0.284179 *10**1,
                0.741609 * 10**1,
                0.891844 * 10**3,
                -0.161309 *10**4,
                0.622106 * 10 **3,
                -0.207588 * 10**3,
                -0.687393 * 10**1,
                0.350716 *10**1 ]
    Koef_n = [1,4,8,9,12,14,0,1,1,3,3,4,5,2,4,0]
    Koef_m = [0,0,0,0,0,0,1,1,2,3,5,5,5,6,6,8]
    

    ## Calculation

    sum = 0
    for i in range(16):
        sum = sum + (Koef_a[i]*((T/T_0)-1)**Koef_m[i])*X_mol_NH3**Koef_n[i]
    
    h = h_0*sum
    
    return h*10**3 # Unit conversion from kJ/kg zu J/kg

def Calc_h_liquid_from_p_X_NH3H2O_Patek( p, X_mol_NH3):
    """
    % ---------------------------------------------------------------------- %
    % Calc_h_from_T_X_NH3H2O_Patek
    % Uses coefficients and formula obtained from Patek 1995
    % ---------------------------------------------------------------------- %
    %{
    Author  : Ludwig Irrgang
    Date    : 25.06.2022
    Copyright information:
    Ludwig Irrgang
    Lehrstuhl für Energiesysteme
    TUM School of Engineering and Design
    Technische Universität München
    Boltzmannstr. 15 
    85748 Garching b. München
    ludwig.irrgang@tum.de
    %}
    % ---------------------------------------------------------------------- %
    % Input:
    %       -   Termperature of Solution T                              [K]
    %       -   Molar Concentration X of NH3 in Solution                [-]
    % Output:
    %       -   Specific enthalpy h of solution                         [J/kg]
    % ---------------------------------------------------------------------- %
    """
    ## Constants
    h_0 = 100                   # [kJ/kg]
    T_0 = 273.16                   #[K]
    T = Calc_T_from_p_X_NH3H2OSol_Patek(p, X_mol_NH3)

    # Table 4
    Koef_a=[  -0.761080 * 10**1,
                0.256905 * 10**2,
                -0.247092 * 10**3,
                0.325952 * 10**3,
                -0.158854 * 10**3,
                0.619084 * 10**2,
                0.114314 * 10**2,
                0.118157 * 10**1,
                0.284179 *10**1,
                0.741609 * 10**1,
                0.891844 * 10**3,
                -0.161309 *10**4,
                0.622106 * 10 **3,
                -0.207588 * 10**3,
                -0.687393 * 10**1,
                0.350716 *10**1 ]
    Koef_n = [1,4,8,9,12,14,0,1,1,3,3,4,5,2,4,0]
    Koef_m = [0,0,0,0,0,0,1,1,2,3,5,5,5,6,6,8]
    

    ## Calculation

    sum = 0
    for i in range(16):
        sum = sum + (Koef_a[i]*((T/T_0)-1)**Koef_m[i])*X_mol_NH3**Koef_n[i]
    
    h = h_0*sum
    
    return h*10**3 # Unit conversion from kJ/kg zu J/kg

def Calc_h_gas_from_T_Y_NH3H2O_Patek(T,Y_mol_NH3):
    """
    % ---------------------------------------------------------------------- %
    % Calc_h_from_T_X_NH3H2O_Patek
    % Uses coefficients and formula obtained from Patek 1995
    % ---------------------------------------------------------------------- %
    %{
    Author  : Ludwig Irrgang
    Date    : 25.06.2022
    Copyright information:
    Ludwig Irrgang
    Lehrstuhl für Energiesysteme
    TUM School of Engineering and Design
    Technische Universität München
    Boltzmannstr. 15 
    85748 Garching b. München
    ludwig.irrgang@tum.de
    %}
    % ---------------------------------------------------------------------- %
    % Input:
    %       -   Termperature of Solution T                              [K]
    %       -   Molar Concentration X of NH3 in Solution                [-]
    % Output:
    %       -   specific enthalpy h of solution                         [J/kg]
    % ---------------------------------------------------------------------- %
    """
    ## Constants
    h_0 = 1000                   # [kJ/kg]
    T_0 = 324                   #[K]
    # Table 5
    Koef_a=[  0.128827 * 10**1,
                0.125247*10**0,
                -0.208748*10**1,
                0.217696*10**1,
                0.235687*10**1,
                -0.886987*10**1,
                0.102635*10**2,
                -0.237440*10**1,
                -0.670514*10**1,
                0.164508*10**2,
                -0.936849*10**1,
                0.842254*10**1,
                0.858807*10**1,
                0.277049*10**1,
                -0.961248*10**0,
                0.988009*10**0,
                0.308482*10**0 ]
    Koef_n = [0,0,0,0,2,2,2,2,3,3,3,4,4,5,6,7,10]
    Koef_m = [0,1,2,3,0,1,2,3,0,1,2,0,1,0,4,2,1]
    

    ## Calculation

    sum = 0
    for i in range(17):
        sum = sum + (Koef_a[i]*(1-(T/T_0))**Koef_m[i])*(1-Y_mol_NH3)**(0.25*Koef_n[i])
    
    h = h_0*sum
    
    return h*10**3 # Unit conversion from kJ/kg zu J/kg

def Calc_h_gas_from_p_X_NH3H2O_Patek(p, X_mol_NH3):
    """
    % ---------------------------------------------------------------------- %
    % Calc_h_from_T_X_NH3H2O_Patek
    % Uses coefficients and formula obtained from Patek 1995
    % ---------------------------------------------------------------------- %
    %{
    Author  : Ludwig Irrgang
    Date    : 25.06.2022
    Copyright information:
    Ludwig Irrgang
    Lehrstuhl für Energiesysteme
    TUM School of Engineering and Design
    Technische Universität München
    Boltzmannstr. 15 
    85748 Garching b. München
    ludwig.irrgang@tum.de
    %}
    % ---------------------------------------------------------------------- %
    % Input:
    %       -   Termperature of Solution T                              [K]
    %       -   Molar Concentration X of NH3 in Solution                [-]
    % Output:
    %       -   specific enthalpy h of solution                            [J/kg]
    % ---------------------------------------------------------------------- %
    """
    ## Constants
    h_0 = 1000                   # [kJ/kg]
    T_0 = 324                    # [K]
    Y_mol_NH3 = Calc_Y_from_p_X_NH3H2OSol_Patek(p,X_mol_NH3)
    T =  Calc_T_from_p_Y_NH3H2OSol_Patek(p, Y_mol_NH3)
    # Table 5
    Koef_a=[  0.128827 * 10**1,
                0.125247*10**0,
                -0.208748*10**1,
                0.217696*10**1,
                0.235687*10**1,
                -0.886987*10**1,
                0.102635*10**2,
                -0.237440*10**1,
                -0.670514*10**1,
                0.164508*10**2,
                -0.936849*10**1,
                0.842254*10**1,
                -0.858807*10**1,
                -0.277049*10**1,
                -0.961248*10**0,
                0.988009*10**0,
                0.308482*10**0 ]
    Koef_n = [0,0,0,0,2,2,2,2,3,3,3,4,4,5,6,7,10]
    Koef_m = [0,1,2,3,0,1,2,3,0,1,2,0,1,0,4,2,1]
    

    ## Calculation

    sum = 0
    for i in range(17):
        sum = sum + (Koef_a[i]*(1-(T/T_0))**Koef_m[i])*(1-Y_mol_NH3)**(0.25*Koef_n[i])
    
    h = h_0*sum
    
    return h*10**3 # Unit conversion from kJ/kg zu J/kg

def Calc_p_from_T_X(T, X, p_start = 400_000):
    '''
    Calculation of:
        p       Pressure of the solution        [Pa]

    Based on:
        X       Molar Concentration X of NH3 in Solution
        T       Temperature of the solution     [K]
    
    Solves the inverse funciton numerical, by using the newton method
    '''
    
    problem = lambda  p,T,X: (Calc_T_from_p_X_NH3H2OSol_Patek(p,X)-T)
    try:
        p_fit = ridder(problem, 100, 4_000_000, args=(T,X))
    except:
        p_fit = bisect(problem, 100, 4_000_000, args=(T,X))

    return p_fit

def NH3inSolution_Calc_X_PT(pressure, Temp):
    """
    % ---------------------------------------------------------------------- %
    % NH3inSolution_Calc_X_PT
    % Since the solution is heated above saturation, some refrigerant must be
    % This function calculates the concentration of the saturated NH3-H2O
    % solution from pressure and temperature
    % Input units: Pressure [Pa], Tempreature [K]
    % ---------------------------------------------------------------------- %
    """
    ## Caluclation
    i = 0
    a = 0.999
    b = 0.001
    try:
        minPressure = Calc_p_from_T_X(Temp,NH3_conversion_w_X("X", b))
        #minPressure = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="TQ", hOut= "P", iUnits=2,iMass=1, iFlag=0, a=Temp, b=0, z=[b, (1-b)]).Output[0] *1000
    except:
        print("Calculation failed in two phase region - only water in solution")
        wNH3 = 0 # only water present in solution
        return wNH3
    
    if minPressure > pressure:
        print('Pressure not in two phase region')
        wNH3 = 0
    else:
        f = lambda w : pressure - Calc_p_from_T_X(Temp,NH3_conversion_w_X("X", w)) #- RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="TQ", hOut= "P", iUnits=2,iMass=1, iFlag=0, a=Temp, b=0, z=[w, (1-w)]).Output[0]*1000 #
        while abs(f(b))> 10: # coresponds to 0.1 kPa
            if f((a+b)/2)<0: # half slitting technique
                a = (a+b)/2
            else:
                b = (a+b)/2
            i = i + 1
            if (i>1000):
                print('Concentration can not be calculated - loop exceeded limits')
        wNH3 = b
    
    return wNH3

def NH3inSolution_Calc_X_PT_REFPROP(pressure, Temp):
    """
    % ---------------------------------------------------------------------- %
    % NH3inSolution_Calc_X_PT
    % Since the solution is heated above saturation, some refrigerant must be
    % This function calculates the concentration of the saturated NH3-H2O
    % solution from pressure and temperature
    % Input units: Pressure [Pa], Tempreature [K]
    % ---------------------------------------------------------------------- %
    """
    ## Caluclation
    i = 0
    a = 0.999
    b = 0.001
    try:
        #minPressure = Calc_p_from_T_X(Temp,NH3_conversion_w_X("X", b))
        minPressure = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="TQ", hOut= "P", iUnits=2,iMass=1, iFlag=0, a=Temp, b=0, z=[b, (1-b)]).Output[0] *1000
    except:
        print("Calculation failed in two phase region - only water in solution")
        wNH3 = 0 # only water present in solution
        return wNH3
    
    if minPressure > pressure:
        print('Pressure not in two phase region')
        wNH3 = 0
    else:
        f = lambda w : pressure - RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="TQ", hOut= "P", iUnits=2,iMass=1, iFlag=0, a=Temp, b=0, z=[w, (1-w)]).Output[0]*1000 #
        while abs(f(b))> 10: # coresponds to 0.1 kPa
            if f((a+b)/2)<0: # half slitting technique
                a = (a+b)/2
            else:
                b = (a+b)/2
            i = i + 1
            if (i>1000):
                print('Concentration can not be calculated - loop exceeded limits')
                sys.exit()
        wNH3 = b
    
    return wNH3

def NH3inSolution_Calc_X_HT(Enthalpy, Temp):
    """
    % ---------------------------------------------------------------------- %
    % NH3inSolution_Calc_X_HT
    % Since the solution is heated above saturation, some refrigerant must be
    % This function calculates the concentration of the saturated NH3-H2O
    % solution from enthalpy and temperature
    % Input units: Enthalpy [J/kg], Tempreature [K]
    % ---------------------------------------------------------------------- %
    """
    ## Caluclation
    i = 0
    a = 0.999
    b = 0.001
    try:
        minEnthalpy = Calc_h_liquid_from_T_X_NH3H2O_Patek(Temp,NH3_conversion_w_X("X", b))
    except:
        print("Calculation failed in two phase region - only water in solution")
        wNH3 = 0 # only water present in solution
    
    if minEnthalpy > Enthalpy:
        print('Enthalpy not in two phase region')
        wNH3 = 0
    else:
        f = lambda w : Enthalpy - Calc_h_liquid_from_T_X_NH3H2O_Patek(Temp,NH3_conversion_w_X("X", w))
        while math.abs(f(b))> 0.001: # coresponds to 1 Pa
            if f((a+b)/2)<0: # half slitting technique
                a = (a+b)/2
            else:
                b = (a+b)/2
            i = i + 1
            if (i>1000):
                print('Concentration can not be calculated - loop exceeded limits')
        wNH3 = b
    
    return wNH3

