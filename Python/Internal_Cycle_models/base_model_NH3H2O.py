import CoolProp.CoolProp as CoolProp
from Fluids import NH3H2O
import numpy
import math
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import os

os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])

RP_Unit = RP.GETENUMdll(iFlag=0,hEnum='MASS SI').iEnum
print(RP_Unit)


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

def base_model_NH3H2O(T, p, h, m, eta, Q, HX, s):
    """
    ## Function Base Model NH3H2O
    # ----------------------------------------------------------------------- #
    #{
    Author  : Ludwig Irrgang
    Date    : 01.09.2022
    Copyright information:
    Ludwig Irrgang
    Lehrstuhl für Energiesysteme
    TUM School of Engineering and Design
    Technische Universität München
    Boltzmannstr. 15 
    85748 Garching b. München
    ludwig.irrgang@tum.de
    #}
    # ----------------------------------------------------------------------- #
   
    if nargin<8||isempty(s),error('Input Argument: Setup s missing')
    if nargin<7||isempty(HX),error('Input Argument: Approach temperature missing')
    if nargin<6||isempty(Q),error('Input Argument: Heat missing')
    if nargin<5||isempty(eta),error('Input Argument: Efficiency missing')
    if nargin<4||isempty(m),error('Input Argument: Mass missing')
    if nargin<3||isempty(h),error('Input Argument: Enthalpie missing')
    if nargin<2||isempty(p),error('Input Argument: Pressure missing')
    if nargin<1||isempty(T),error('Input Argument: Temperature missing')
    # ----------------------------------------------------------------------- #
    ## Input/ Output
    # ----------------------------------------------------------------------- #
    # Input:
    #{
    Initialize structs with:
    - T.evap
    - T.sol_abs_out
    - T.sol_des_out
    - T.cond
    - T.ext_cond_in (refrigerant is subcooled)
    - Q.dec
    - eta.pump
    - HX.T_PP_SHEX
    - HX.T_PP_RHEX
    - HX.T_PP_cond
    - HX.dT_ref_des
    #}
    # Output:
    #{
    1. Temperature struct                       [K]
    2. Pressure struct                          [Pa]
    3. Specific enthalpy struct                 [J/kg]
    4. Mass flow rate struct                    [kg/s]
    5. Mass fraction struct                     [kg/s]
    6. Efficiency struct                        [-]
    7. Heat flow struct                         [J]
    8. Post process struct
    9. Setup/ Entropy struct
    #}
    # ----------------------------------------------------------------------- #
    ## Absorption System
    # ----------------------------------------------------------------------- #
    # Components:
    #{
    - Desorber
    - Condenser
    - Evaporator
    - Absorber
    - Pump
    - 2 x Throttle Valves
    - Solution Heat Exchanger
    - Refrigerant Heat Exchanger
    - Working fluid: Ammonia Water Solution
    - Refrigerant: Ammonia
    #}
    # ----------------------------------------------------------------------- #
    ## Assumptions
    # ----------------------------------------------------------------------- #
    #{
    - Temperature of refrigerant leaving desorber is 5K below desorber temp.
    - All components of the system operate in steady state
    - Solution leaving generator and absorber is saturated
    - Refrigerant leaving condenser and evaporator is saturated
    - Pressure drops in the system components are negelcted
    #} 
    """

    ## Initialise requierd variables
    class var:
        i = 0 
    w = var()
    x = var()
    v= var()
    PP= var()
    rho = var()
    cp = var()
    # -----------------------------------------------------------------------  #
    # # Definition of constants
    # ------------------------Necessary Constants----------------------------  #
    M_NH3 = 0.017031                 #[kg/mol]
    M_H2O = 0.018015268              #[kg/mol]
    # ----------------------------------------------------------------------- #
    ## Calculation
    # ---------------------------CALCULATION--------------------------------- #
    ## Calculate pressures
    # Low pressure
    p.evap = CoolProp.PropsSI('P','T',T.evap,'Q',1,'AMMONIA')
    p.cond = CoolProp.PropsSI('P','T',T.cond,'Q',0,'AMMONIA')
    ## Refrigerant line
    # Desorber
    T.ref_des_out = T.sol_des_out - HX.dT_ref_des
    h.ref_des_out = CoolProp.PropsSI('H','T',T.ref_des_out,'P',p.cond,'AMMONIA')
    # Condenser
    T.ref_cond_in = T.ref_des_out
    h.ref_cond_in = h.ref_des_out
    if (T.ext_cond_in+HX.T_PP_cond<T.cond):
        T.ref_cond_out = T.ext_cond_in + HX.T_PP_cond # Subcooling as low as possible
        h.ref_cond_out = CoolProp.PropsSI('H','T',T.ref_cond_out,'P',p.cond,'AMMONIA')
    else:
        T.ref_cond_out = T.cond
        h.ref_cond_out = CoolProp.PropsSI('H','T',T.ref_cond_out,'Q',0,'AMMONIA')
    
    # Evaporator
    T.ref_evap_out = T.evap
    h.ref_evap_out = CoolProp.PropsSI('H','P',p.evap,'Q',1,'AMMONIA')
    # Subcooler (heat capacity of steam lower than liquid)
    if(T.ref_cond_out - HX.T_PP_RHEX > T.ref_evap_out):
        T.ref_abs_in = T.ref_cond_out - HX.T_PP_RHEX
        h.ref_abs_in = CoolProp.PropsSI('H','T',T.ref_abs_in,'P',p.evap,'AMMONIA')
    else:
        T.ref_abs_in = T.ref_evap_out
        h.ref_abs_in = h.ref_evap_out
    
    h.ref_valve_in = h.ref_cond_out - (h.ref_abs_in-h.ref_evap_out)
    T.ref_valve_in = CoolProp.PropsSI('T','H',h.ref_valve_in,'P',p.cond,'AMMONIA')
    # Throttle
    h.ref_evap_in = h.ref_valve_in # Isenthalpic throttle
    T.ref_evap_in = CoolProp.PropsSI('T','H',h.ref_evap_in,'P',p.evap,'AMMONIA')
    #-------------------------------------------------------------------------#
    ## Rich solution (High ref. concentration)
    # Absorber
    w.NH3_rich = NH3H2O.NH3inSolution_Calc_X_PT(p.evap,T.sol_abs_out)
    w.NH3_rich_2 = w.NH3_rich
    w.NH3_rich = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="PT", hOut= "XMASS", iUnits=RP_Unit,iMass=1, iFlag=0, a=p.evap*10**(-6), b=T.sol_abs_out, z=[w.NH3_rich, (1-w.NH3_rich)]).Output[0]
    x.NH3_rich = NH3_conversion_w_X("X",w.NH3_rich)
    x.H2O_rich = 1-x.NH3_rich
    h.sol_abs_out = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="TQ", hOut= "H", iUnits=RP_Unit,iMass=1, iFlag=0, a=T.sol_abs_out, b=0, z=[w.NH3_rich, (1-w.NH3_rich)]).Output[0]*1000
    # Pump
    s.sol_abs_out = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="TQ", hOut= "S", iUnits=RP_Unit,iMass=1, iFlag=0, a=T.sol_abs_out, b=0, z=[w.NH3_rich, (1-w.NH3_rich)]).Output[0]*1000
    s.sol_pump_out = s.sol_abs_out
    h.sol_pump_isentropic = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="PS", hOut= "H", iUnits=RP_Unit,iMass=1, iFlag=0, a=p.cond*10**(-6), b=s.sol_pump_out*10**(-3), z=[w.NH3_rich, (1-w.NH3_rich)]).Output[0]*1000
    h.sol_pump_out = h.sol_abs_out - (h.sol_abs_out-h.sol_pump_isentropic)/eta.pump
    T.sol_pump_out = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="PH", hOut= "T", iUnits=RP_Unit, iMass=1, iFlag=0, a=p.cond*10**(-6), b=h.sol_pump_out*10**(-3), z=[w.NH3_rich, (1-w.NH3_rich)]).Output[0]
    cp.sol_pump_out = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="PH", hOut= "CP", iUnits=RP_Unit, iMass=1, iFlag=0, a=p.cond*10**(-6), b=h.sol_pump_out*10**(-3), z=[w.NH3_rich, (1-w.NH3_rich)]).Output[0]*1000
    #-------------------------------------------------------------------------#
    ## Poor solution (Low ref. concentration)
    # Desorber
    w.NH3_poor = NH3H2O.NH3inSolution_Calc_X_PT(p.cond,T.sol_des_out)
    w.NH3_poor = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="PT", hOut= "XMASS", iUnits=RP_Unit,iMass=1, iFlag=0, a=p.cond*10**(-6), b=T.sol_des_out, z=[w.NH3_poor, (1-w.NH3_poor)]).Output[0]
    if(w.NH3_poor>0):
        h.sol_des_out = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="TQ", hOut= "H", iUnits=RP_Unit,iMass=1, iFlag=0, a=T.sol_des_out, b=0, z=[w.NH3_poor, (1-w.NH3_poor)]).Output[0]*1000
    else:
        h.sol_des_out = CoolProp.PropsSI('H', 'T', T.sol_des_out, 'P', p.cond, 'WATER')
    
    # SHEX
    if(T.sol_pump_out + HX.T_PP_SHEX < T.sol_des_out):
        T.sol_valve_in = T.sol_pump_out + HX.T_PP_SHEX
        T.sol_valve_in = T.sol_des_out
    else:
        T.sol_valve_in = T.sol_des_out

    h.sol_valve_in = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="PT", hOut= "H", iUnits=RP_Unit,iMass=1, iFlag=0, a=p.cond*10**(-6), b=T.sol_valve_in, z=[w.NH3_poor, (1-w.NH3_poor)]).Output[0]*1000
    #-------------------------------------------------------------------------#
    ## Energy balances and mass conservation
    # Solve system of linear equations (energy conservation, mass conservation)
    # Solution vector: m_rich, m_ref, m_poor
    match s.requirement:
        case 'Q_des':
            A = numpy.array([[   h.sol_pump_out,     -h.ref_des_out,     -h.sol_vavle_in],
                    [w.NH3_rich,         -1,                 -w.NH3_poor],
                    [1,                  -1,                 -1   ]])
            b = numpy.array([   -Q.dec,     0,      0   ])
        case 'Q_evap':
            A = numpy.array([[   0,          (h.ref_evap_out-h.ref_evap_in), 0 ],
                    [w.NH3_rich, -1,                             -w.NH3_poor],
                    [1,          -1,                             -1   ]])
            b = numpy.array([   Q.dec  ,   0,      0   ])
        case _:
            print('AKM requirement is not defined properly. Use Q_des or Q_evap')
             # error('AKM requirement is not defined properly. Use Q_des or Q_evap')
    try:
        y = numpy.linalg.solve(A,b)
        m.sol_rich = y[0]
        m.ref = y[1]
        m.sol_poor = y[2]
        print("linalg.solve")
    except:
        y = numpy.linalg.lstsq(A,b)[0]
        m.sol_rich = y[0]
        m.ref = y[1]
        m.sol_poor = y[2]
        print("linalg.lstsq")

    #-------------------------------------------------------------------------#
    ## Check
    # Refrigerant concentrations
    if (w.NH3_rich < w.NH3_poor):
        print("w_NH3_rich  < w_NH3_poor ")
        #error("w_NH3_rich < w_NH3_poor")
    
    if (w.NH3_rich < 0):
        print("w_NH3_rich < 0")
        #error("w_NH3_rich < 0")
    
    if (w.NH3_poor < 0):
        print("w_NH3_poor < 0")
        #error("w_NH3_poor < 0")
    
    if (w.NH3_rich - w.NH3_poor < 0.001):
        print("w_NH3_rich - w_NH3_poor < 0.001")
        #error("w_NH3_rich - w_NH3_poor < 0.001")
    
    # Mass flows
    if (m.ref<0 or m.sol_poor<0 or m.sol_rich<0):
        print("mass flow is negativ")
       # error("mass flow is negativ")
    
    # SHEX
    h.sol_des_in = (m.sol_poor*h.sol_des_out + m.sol_rich*h.sol_pump_out - m.sol_poor*h.sol_valve_in) / m.sol_rich
    try:
        T.sol_des_in = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="PH", hOut= "T", iUnits=RP_Unit,iMass=1, iFlag=0, a=p.cond*10**(-6), b=h.sol_des_in*10**(-3), z=[w.NH3_rich, (1-w.NH3_rich)]).Output[0]
    except:
        # Refprop fails for specific enthalpy values
        T.sol_des_in = T.sol_pump_out + (h.sol_des_in - h.sol_pump_out)/cp.sol_pump_out

    # Valve
    h.sol_abs_in = h.sol_valve_in
    try: 
        T.sol_abs_in = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="PH", hOut= "T", iUnits=RP_Unit,iMass=1, iFlag=0, a=p.evap*10**(-6), b=h.sol_abs_in*10**(-3), z=[w.NH3_poor, (1-w.NH3_poor)]).Output[0]
    except:
        # Refprop fails for sepcific enthalpy values
        T.sol_abs_in = T.sol_valve_in
    #-------------------------------------------------------------------------#
    ## Post processing
    # Fluxes over system boundary
    Q.cond = m.ref*(h.ref_cond_out-h.ref_cond_in)
    Q.evap = m.ref*(h.ref_evap_out-h.ref_evap_in)
    Q.abs =  m.sol_rich*h.sol_abs_out - m.ref*h.ref_abs_in - m.sol_poor*h.sol_abs_in
    PP.W_pump = m.sol_rich*(h.sol_pump_out-h.sol_abs_out)
    Q.des = m.ref*h.ref_des_out + m.sol_poor*h.sol_des_out - m.sol_rich*h.sol_des_in
    # Heat Exchanger
    Q.SHEX = m.sol_poor*(h.sol_des_out-h.sol_valve_in)
    Q.RHEX = m.ref*(h.ref_abs_in-h.ref_evap_out)
    h.RHEXideal = CoolProp.PropsSI('H','T',T.ref_cond_out,'P',p.evap,'AMMONIA')
    eta.RHEX = (h.ref_abs_in-h.ref_evap_out)/(h.RHEXideal-h.ref_evap_out)
    eta.SHEX = (h.sol_des_in-h.sol_pump_out)/(h.sol_des_out-h.sol_pump_out)
    # COP
    PP.COP = Q.evap/(Q.des + PP.W_pump)
    # Circulation
    PP.f = m.sol_rich/m.ref
    # Energy balance
    PP.energyBalance = Q.des + Q.evap + PP.W_pump + Q.cond + Q.abs
    # Mass balance
    PP.massBalance = m.ref + m.sol_poor - m.sol_rich
    
    #-------------------------------------------------------------------------#
    ## Check
    # Energy and mass balance
    if (abs(PP.energyBalance) > 1):
        print("Energy is not conserved")
      #  error("Energy is not conserved")
    
    if (abs(PP.massBalance) > 1):
        print("Mass is not conserved")
       # error("Mass is not conserved")
    
    #-------------------------------------------------------------------------#
    return T, p, h, m, w, eta, Q, PP, s
    