import CoolProp.CoolProp as CoolProp
from Fluids import NH3H2O
import numpy
import math

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
    T.ref_des_out = T.sol_des_out - 5
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
    else:
        T.ref_abs_in = T.ref_evap_out
    
    h.ref_abs_in = CoolProp.PropsSI('H','T',T.ref_abs_in,'P',p.evap,'AMMONIA')
    h.ref_valve_in = h.ref_cond_out - (h.ref_abs_in-h.ref_evap_out)
    T.ref_valve_in = CoolProp.PropsSI('T','H',h.ref_valve_in,'P',p.cond,'AMMONIA')
    # Throttle
    h.ref_evap_in = h.ref_valve_in # Isenthalpic throttle
    T.ref_evap_in = CoolProp.PropsSI('T','H',h.ref_evap_in,'P',p.evap,'AMMONIA')
    #-------------------------------------------------------------------------#
    ## Rich solution (High ref. concentration)
    # Absorber
    w.NH3_rich = NH3H2O.NH3inSolution_Calc_X_PT(p.evap,T.sol_abs_out)
    x.NH3_rich = NH3_conversion_w_X("X",w.NH3_rich)
    x.H2O_rich = 1-x.NH3_rich
   
    h.sol_abs_out = NH3H2O.Calc_h_liquid_from_T_X_NH3H2O_Patek(T.sol_abs_out, NH3_conversion_w_X("X", w.NH3_rich))
    # Pump
    rho.NH3_molar_sol_abs_out = CoolProp.PropsSI("DMOLAR","T", T.sol_abs_out, "P", p.cond,"AMMONIA") # [mol/m3]
    rho.H2O_molar_sol_abs_out = CoolProp.PropsSI("DMOLAR","T", T.sol_abs_out, "P", p.cond,"WATER")   # [mol/m3]
    v.mol_sol_abs_out = (1/rho.NH3_molar_sol_abs_out)*x.NH3_rich + (1/rho.H2O_molar_sol_abs_out)*x.H2O_rich # [m3/mol]
    M_mix = x.NH3_rich*M_NH3 + x.H2O_rich*M_H2O # [kg/mol]
    v.sol_abs_out = v.mol_sol_abs_out/M_mix     # [m3/kg]    
    w.pump = v.sol_abs_out*(p.cond - p.evap)
    h.sol_pump_out = h.sol_abs_out + w.pump/eta.pump
    T.sol_pump_out = T.sol_abs_out # copied from LiBr, assuming isothermal

    # rho.sol_abs_out_mol = NH3H2O.Calc_rho_from_T_X
    # rho.sol_abs_out = rho.sol_abs_out_mol * (NH3_conversion_w_X("X",w.NH3_rich)*M_NH3 +(1-NH3_conversion_w_X("X",w.NH3_rich))*M_H2O)
    # v.sol_abs_out = 1/rho.sol_abs_out 
    # w.pump = v.sol_abs_out*(p.cond-p.evap)
    # h.sol_pump_out = h.sol_abs_out + w.pump/eta.pump 
    # s.sol_abs_out = refpropm('S','T',T.sol_abs_out,'Q',0,'AMMONIA','WATER',[w.NH3_rich (1-w.NH3_rich)])
    # s.sol_pump_out = s.sol_abs_out
    #h.sol_pump_isentropic = refpropm('H','P',p.cond/1000,'S',s.sol_pump_out,'AMMONIA','WATER',[w.NH3_rich (1-w.NH3_rich)])
    #h.sol_pump_out = h.sol_abs_out - (h.sol_abs_out-h.sol_pump_isentropic)/eta.pump
    #T.sol_pump_out = refpropm('T','P',p.cond/1000,'H',h.sol_pump_out,'AMMONIA','WATER',[w.NH3_rich (1-w.NH3_rich)])
    #-------------------------------------------------------------------------#
    ## Poor solution (Low ref. concentration)
    # Desorber
    w.NH3_poor = NH3H2O.NH3inSolution_Calc_X_PT(p.cond,T.sol_des_out)
    if(w.NH3_poor>0):
        h.sol_des_out = NH3H2O.Calc_h_liquid_from_T_X_NH3H2O_Patek(T.sol_des_out, NH3_conversion_w_X("X", w.NH3_poor))
    else:
        h.sol_des_out = CoolProp.PropsSI('H', 'T', T.sol_des_out, 'P', p.cond, 'WATER')
    
    # SHEX
    if(T.sol_pump_out + HX.T_PP_SHEX < T.sol_des_out):
        T.sol_valve_in = T.sol_pump_out + HX.T_PP_SHEX
    else:
        T.sol_valve_in = T.sol_des_out
    
 
    h.sol_valve_in = NH3H2O.Calc_h_liquid_from_p_X_NH3H2O_Patek(p.cond, NH3_conversion_w_X("X", w.NH3_poor)) # h_liquid or h_gas?
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
    #hier evtl mit der Determinante mal vorher prüfen.
    d= numpy.linalg.det(A)
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
    
    # SHEX
    h.sol_des_in = (m.sol_poor*h.sol_des_out + m.sol_rich*h.sol_pump_out - m.sol_poor*h.sol_valve_in) / m.sol_rich
    T.sat_rich_SHEX = NH3H2O.Calc_T_from_p_X_NH3H2OSol_Patek(p.cond, NH3_conversion_w_X("X",w.NH3_rich))
    h.sat_rich_SHEX = NH3H2O.Calc_h_liquid_from_p_X_NH3H2O_Patek(p.cond, NH3_conversion_w_X("X",w.NH3_rich))
    if (h.sat_rich_SHEX>h.sol_des_in): # Saturation is not reached in SHEX - no evaporation
        T.sol_des_in = refpropm('T','P',p.cond/1000,'H',h.sol_des_in,'AMMONIA','WATER',[w.NH3_rich (1-w.NH3_rich)])
    else:
        T.sol_des_in = NH3H2O.NH3inSolution_Calc_state_SHEX_exit(T.sat_rich_SHEX, p.cond, h.sol_des_in, m.sol_rich, w.NH3_rich, 0.1)
    
    # Valve
    h.sol_abs_in = h.sol_valve_in
    T.sat_poor_valve = NH3H2O.Calc_T_from_p_X_NH3H2OSol_Patek(p.evap, NH3_conversion_w_X("X",w.NH3_poor))
    h.sat_poor_valve = NH3H2O.Calc_h_liquid_from_p_X_NH3H2O_Patek(p.evap, NH3_conversion_w_X("X",w.NH3_poor))
    if (h.sat_poor_valve>h.sol_abs_in): # Saturation is not reached after valve - no evaporation
        T.sol_abs_in = T.sol_valve_in
    else:
        T.sol_abs_in = NH3H2O.NH3inSolution_Calc_state_valve_exit(T.sat_poor_valve, p.evap, h.sol_abs_in, m.sol_poor, 1-w.NH3_poor, 0.1)
    
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
    
    # Energy and mass balance
    if (math.abs(PP.energyBalance) > 1):
        print("Energy is not conserved")
      #  error("Energy is not conserved")
    
    if (math.abs(PP.massBalance) > 1):
        print("Mass is not conserved")
       # error("Mass is not conserved")
    
    #-------------------------------------------------------------------------#
    return T, p, h, m, w, eta, Q, PP, s
    