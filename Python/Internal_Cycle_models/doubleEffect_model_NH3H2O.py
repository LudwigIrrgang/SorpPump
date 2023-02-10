import CoolProp.CoolProp as CoolProp
import math
import numpy
from Fluids import LiBrSol, NH3H2O, H2O

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

def doubleEffect_model_NH3H2O(T, p, h, m, eta, Q, HX, s):
    """
    ## Function Double Effect Model NH3H2O (Parallel Flow Configuration)
    # ----------------------------------------------------------------------- #
    #{
    Author  : Ludwig Irrgang
    Date    : 02.09.2022
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
    - T.cond_int
    - Q.dec
    - eta.pump
    - HX.T_PP_SHEX
    - HX.T_PP_SHEXI
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
    9. Setup struct
    #}
    # ----------------------------------------------------------------------- #
    ## Absorption System
    # ----------------------------------------------------------------------- #
    # Components:
    #{
    - 2 x Desorber
    - 2 x Condenser
    - 2 x Pump
    - 3 x Throttle Valves
    - Evaporator
    - Absorber
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
    - Solution leaving desorber and absorber is saturated
    - Refrigerant leaving condenser and evaporator is saturated
    - Pressure drops in the system components are negelcted
    - Heat capacity of solution assumed to be constant in SHEX
    #} 
    """
    # -----------------------------------------------------------------------  #
    # # Initialisation of variable classes for storing calculatioin results
    # -----------------------------------------------------------------------  #
    class var:
        i = 0 
    
    w= var()
    PP = var()
    # ----------------------------------------------------------------------- #
    ## CalculationCoolProp
    # ---------------------------CALCULATION--------------------------------- #
    ## Calculate pressures
    p.evap = CoolProp.PropsSI('P','T',T.evap,'Q',1,'Ammonia')
    p.cond = CoolProp.PropsSI('P','T',T.cond,'Q',0,'Ammonia')
    p.cond_int = CoolProp.PropsSI('P','T',T.cond_int,'Q',0,'Ammonia')
    # ----------------------------------------------------------------------- #
    ## Find minimal possible internal condensation temperature
    w.NH3_rich = NH3H2O.NH3inSolution_Calc_X_PT(p.evap/1000,T.sol_abs_out)
    w.NH3_poor = NH3H2O.NH3inSolution_Calc_X_PT(p.cond/1000,T.sol_des_out)
    if(w.NH3_rich<w.NH3_poor):
        w.NH3_poor = w.NH3_rich-0.05 # Ausgasungsbreite auf 5 # eingestellt
        T.sol_des_out = NH3H2O.Calc_T_from_p_X_NH3H2OSol_Patek(p.cond*10**(-6), NH3_conversion_w_X("X", w.NH3_poor))
        T.cond_int = T.sol_des_out + HX.T_PP_cond_int
        p.cond_int = CoolProp.PropsSI('P','T',T.cond_int,'Q',0,'Ammonia')

    ## Refrigerant line
    # Desorber
    T.ref_des_outI = T.sol_des_outI - 5
    h.ref_des_outI = CoolProp.PropsSI('H','T',T.ref_des_outI,'P',p.cond_int,'Ammonia')
    T.ref_des_out = T.sol_des_out - 5
    h.ref_des_out = refpropm('H','T',T.ref_des_out,'P',p.cond/1000,'Ammonia') # Kann nicht mit CoolProp gerechnet werden
    # Condenser
    T.ref_cond_inI = T.ref_des_outI
    h.ref_cond_inI = h.ref_des_outI
    T.ref_cond_outI = T.cond_int
    h.ref_cond_outI = CoolProp.PropsSI('H','P',p.cond_int,'Q',0,'Ammonia')
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
    h.ref_evap_out = CoolProp.PropsSI('H','P',p.evap,'Q',1,'Ammonia')
    # Subcooler
    if(T.ref_cond_out - HX.T_PP_RHEX > T.ref_evap_out):
        T.ref_abs_in = T.ref_cond_out - HX.T_PP_RHEX
    else:
        T.ref_abs_in = T.ref_evap_out

    h.ref_abs_in = CoolProp.PropsSI('H','T',T.ref_abs_in,'P',p.evap,'Ammonia') 
    h.ref_valve_in = h.ref_cond_out - (h.ref_abs_in-h.ref_evap_out)
    T.ref_valve_in = CoolProp.PropsSI('T','H',h.ref_valve_in,'P',p.cond,'Ammonia')
    # Expansion valve
    h.ref_valve_outI = h.ref_cond_outI # Isenthalpic thottle
    T.ref_valve_outI = CoolProp.PropsSI('T','H',h.ref_valve_outI,'P',p.cond,'Ammonia')
    h.ref_evap_in = h.ref_valve_in # Isenthalpic thottle
    T.ref_evap_in = CoolProp.PropsSI('T','H',h.ref_evap_in,'P',p.evap,'Ammonia')
    # ----------------------------------------------------------------------- #
    ## Rich solution (High ref. concentration)
    # Absorber
    w.NH3_rich = NH3H2O.NH3inSolution_Calc_X_PT(p.evap/1000,T.sol_abs_out)
    h.sol_abs_out = refpropm('H','T',T.sol_abs_out,'P',p.evap/1000,'AMMONIA','WATER',[w.NH3_rich (1-w.NH3_rich)])
    # Pump
    # Low pressure
    # s.sol_abs_out = refpropm('S','T',T.sol_abs_out,'P',p.evap/1000,'AMMONIA','WATER',[w.NH3_rich (1-w.NH3_rich)])
    # s.sol_pump_out = s.sol_abs_out
    # h.sol_pump_isentropic = refpropm('H','P',p.cond/1000,'S',s.sol_pump_out,'AMMONIA','WATER',[w.NH3_rich (1-w.NH3_rich)])
    # h.sol_pump_out = h.sol_abs_out - (h.sol_abs_out-h.sol_pump_isentropic)/eta.pump
    # T.sol_pump_out = refpropm('T','P',p.cond/1000,'H',h.sol_pump_out,'AMMONIA','WATER',[w.NH3_rich (1-w.NH3_rich)])
    # # High pressure
    # h.sol_pump_isentropicI = refpropm('H','P',p.cond_int/1000,'S',s.sol_pump_out,'AMMONIA','WATER',[w.NH3_rich (1-w.NH3_rich)])
    # h.sol_pump_outI = h.sol_abs_out - (h.sol_abs_out-h.sol_pump_isentropicI)/eta.pump
    # T.sol_pump_outI = refpropm('T','P',p.cond_int/1000,'H',h.sol_pump_outI,'AMMONIA','WATER',[w.NH3_rich (1-w.NH3_rich)])
    # ----------------------------------------------------------------------- #
    ## Poor solution (Low ref. concentration)
    # Low pressure
    # Desorber
    w.NH3_poor = NH3H2O.NH3inSolution_Calc_X_PT(p.cond/1000,T.sol_des_out)
    if(w.NH3_poor>0):
        h.sol_des_out = NH3H2O.Calc_h_liquid_from_T_X_NH3H2O_Patek(T.sol_des_out, NH3_conversion_w_X("X", w.NH3_poor))
    else:
        h.sol_des_out = H2O.Calc_h_from_T_p_H2O_Patek(T.sol_des_out, p.cond*10**(-6))

    # SHEX
    if(T.sol_pump_out + HX.T_PP_SHEX < T.sol_des_out):
        T.sol_valve_in = T.sol_pump_out + HX.T_PP_SHEX
    else:
        T.sol_valve_in = T.sol_des_out

    h.sol_valve_in = refpropm('H','T',T.sol_valve_in,'P',p.cond/1000,'AMMONIA','WATER',[w.NH3_poor (1-w.NH3_poor)])
    # High pressure
    # Desober
    w.NH3_poorI = NH3H2O.NH3inSolution_Calc_X_PT(p.cond_int/1000,T.sol_des_outI)
    if(w.NH3_poorI>0):
        h.sol_des_outI = NH3H2O.Calc_h_liquid_from_T_X_NH3H2O_Patek(T.sol_des_outI, NH3_conversion_w_X("X", w.NH3_poorI))
    else:
        h.sol_des_outI = H2O.Calc_h_from_T_p_H2O_Patek(T.sol_des_outI, p.cond_int*10**(-6))

    # SHEXI
    if(T.sol_pump_outI + HX.T_PP_SHEXI < T.sol_des_outI):
        T.sol_valve_inI = T.sol_pump_outI + HX.T_PP_SHEXI
    else:
        T.sol_valve_inI = T.sol_des_outI

    h.sol_valve_inI = refpropm('H','T',T.sol_valve_inI,'P',p.cond_int/1000,'AMMONIA','WATER',[w.NH3_poorI (1-w.NH3_poorI)])
    # ----------------------------------------------------------------------- #
    ## Energy balances and mass conservation
    # Solve system of linear equations (energy conservation, mass conservation)
    # Solution vector: Q.condI, Q.des_mid, m.ref, m.sol_rich, m.sol_poor, m.refI, m.sol_richI, m.sol_poorI
    match s.requirement:
        case 'Q_des':
            A = numpy.array([[   0,          0,  0,                  0,              0,              -h.ref_des_out,                     h.sol_pump_outI,    -h.sol_valve_inI]
                    [0,          0,  0,                  0,              0,              -1,                                 w.NH3_rich,         -w.NH3_poorI]
                    [0,          0,  0,                  0,              0,              -1,                                 1,                  -1]
                    [1,          0,  0,                  0,              0,              (h.ref_cond_inI-h.ref_cond_outI),   0,                  0]
                    [1,          1,  0,                  0,              0,              0,                                  0,                  0]
                    [0,          1,  -h.ref_des_out,     h.sol_pump_out, -h.sol_valve_in,0,                                  0,                  0]
                    [0,          0,  -1,                 w.NH3_rich,     -w.NH3_poor,    0,                                  0,                  0]
                    [0,          0,  -1,                 1,              -1,             0,                                  0,                  0] ])
            b = numpy.array([   -Q.dec,     0,  0,                  0,              0,              0,                                  0,                  0   ])
        case 'Q_evap':
            A = numpy.array([  [ 0,          0,  (h.ref_evap_out-h.ref_evap_in),     0,              0,                  (h.ref_evap_out-h.ref_evap_in),     0,              0]
                    [0,          0,  0,                                  0,              0,                  -1,                                 w.NH3_rich,    -w.NH3_poorI]
                    [0,          0,  0,                                  0,              0,                  -1,                                 1,              -1]
                    [1,          0,  0,                                  0,              0,                  (h.ref_cond_inI-h.ref_cond_outI),   0,              0]
                    [1,          1,  0,                                  0,              0,                  0,                                  0,              0]
                    [0,          1,  -h.ref_des_out,                     h.sol_pump_out, -h.sol_valve_in,    0,                                  0,              0]
                    [0,          0,  -1,                                 w.NH3_rich,     -w.NH3_poor,        0,                                  0,              0]
                    [0,          0,  -1,                                 1,              -1,                 0,                                  0,              0]   ])
            b = numpy.array([   Q.dec,      0,  0  ,                                0  ,            0    ,              0     ,                             0 ,             0   ])
        case _:
            print('AKM decider is not defined properly. Use Q_des or Q_evap')
            #error('AKM decider is not defined properly. Use Q_des or Q_evap')

    y = numpy.linalg.solve(A,b)
    # High pressure side
    Q.condI = y(1)
    m.refI = y(6)
    m.sol_richI = y(7)
    m.sol_poorI = y(8)
    # Low pressure side
    Q.des_mid = y(2)
    m.ref = y(3)
    m.sol_rich = y(4)
    m.sol_poor = y(5)
    ## Rich solution after SHEX
    # Low pressure
    h.sol_des_in = (m.sol_poor*h.sol_des_out + m.sol_rich*h.sol_pump_out - m.sol_poor*h.sol_valve_in) / m.sol_rich
    T.sat_rich_SHEX = NH3H2O.Calc_T_from_p_X_NH3H2OSol_Patek(p.cond*10**(-6),NH3_conversion_w_X("X", w.NH3_rich) )
    h.sat_rich_SHEX = refpropm('H','T',T.sat_rich_SHEX,'P',p.cond/1000,'AMMONIA','WATER',[w.NH3_rich (1-w.NH3_rich)])
    if (h.sat_rich_SHEX>h.sol_des_in): # Saturation is not reached in SHEX - no evaporation
        T.sol_des_in = refpropm('T','P',p.cond/1000,'H',h.sol_des_in,'AMMONIA','WATER',[w.NH3_rich (1-w.NH3_rich)])
    else:
        T.sol_des_in = NH3H2O.NH3inSolution_Calc_state_SHEX_exit(T.sat_rich_SHEX, p.cond, h.sol_des_in, m.sol_rich, w.NH3_rich, 0.1)

    # High pressure
    h.sol_des_inI = (m.sol_poorI*h.sol_des_outI + m.sol_richI*h.sol_pump_outI - m.sol_poorI*h.sol_valve_inI) / m.sol_richI
    T.sat_rich_SHEXI = NH3H2O.Calc_T_from_p_X_NH3H2OSol_Patek(p.cond_int*10**(-6), NH3_conversion_w_X("X", w.NH3_rich))
    h.sat_rich_SHEXI = refpropm('H','T',T.sat_rich_SHEXI,'P',p.cond_int/1000,'AMMONIA','WATER',[w.NH3_rich (1-w.NH3_rich)])
    if (h.sat_rich_SHEXI>h.sol_des_inI): # Saturation is not reached in SHEX - no evaporation
        T.sol_des_inI = refpropm('T','P',p.cond_int/1000,'H',h.sol_des_inI,'AMMONIA','WATER',[w.NH3_rich (1-w.NH3_rich)])
    else:
        T.sol_des_inI = NH3H2O.NH3inSolution_Calc_state_SHEX_exit(T.sat_rich_SHEXI, p.cond_int, h.sol_des_inI, m.sol_richI, w.NH3_rich, 0.1)

    ## Poor solution after valve
    # Low pressure
    h.sol_abs_in = h.sol_valve_in
    T.sat_poor_valve = NH3H2O.Calc_T_from_p_X_NH3H2OSol_Patek(p.evap*10**(-6), NH3_conversion_w_X("X", w.NH3_poor))
    h.sat_poor_valve = refpropm('H','T',T.sat_poor_valve,'P',p.evap/1000,'AMMONIA','WATER',[w.NH3_poor (1-w.NH3_poor)])
    if (h.sat_poor_valve>h.sol_abs_in): # Saturation is not reached after valve - no evaporation
        T.sol_abs_in = T.sol_valve_in
    else:
        T.sol_abs_in = NH3H2O.NH3inSolution_Calc_state_valve_exit(T.sat_poor_valve, p.evap, h.sol_abs_in, m.sol_poor, w.NH3_poor, 0.1)

    # High pressure
    h.sol_abs_inI = h.sol_valve_inI
    T.sat_poor_valveI = NH3H2O.Calc_T_from_p_X_NH3H2OSol_Patek(p.evap*10**(-6), NH3_conversion_w_X("X", w.NH3_poorI))
    h.sat_poor_valveI = refpropm('H','T',T.sat_poor_valveI,'P',p.evap/1000,'AMMONIA','WATER',[w.NH3_poorI (1-w.NH3_poorI)])
    if (h.sat_poor_valveI>h.sol_abs_inI): # Saturation is not reached after valve - no evaporation
        T.sol_abs_inI = T.sol_valve_inI
    else:
        T.sol_abs_inI = NH3H2O.NH3inSolution_Calc_state_valve_exit(T.sat_poor_valveI, p.evap, h.sol_abs_inI, m.sol_poorI, w.NH3_poorI, 0.1)

    # ----------------------------------------------------------------------- #
    ## Post processing
    # Calculate fluxes over system boundary
    Q.cond = h.ref_cond_out*(m.ref+m.refI) - h.ref_cond_in*m.ref - h.ref_valve_outI*m.refI 
    Q.evap = (m.ref+m.refI)*(h.ref_evap_out-h.ref_evap_in)
    Q.abs = (m.sol_rich+m.sol_richI)*h.sol_abs_out - (m.ref+m.refI)*h.ref_abs_in - m.sol_poor*h.sol_abs_in - m.sol_poorI*h.sol_abs_inI
    Q.des =  m.refI*h.ref_des_outI + m.sol_poorI*h.sol_valve_inI - m.sol_richI*h.sol_pump_outI
    PP.W_pump = m.sol_rich*(h.sol_pump_out-h.sol_abs_out)
    PP.W_pumpI = m.sol_richI*(h.sol_pump_outI-h.sol_abs_out)
    # Heat Exchanger
    Q.SHEX = m.sol_poor*(h.sol_des_out-h.sol_valve_in)
    Q.SHEXI = m.sol_poorI*(h.sol_des_outI-h.sol_valve_inI)
    Q.RHEX = m.ref*(h.ref_abs_in-h.ref_evap_out)
    h.RHEXideal = CoolProp.PropsSI('H','T',T.ref_cond_out,'P',p.evap,'AMMONIA')
    eta.RHEX = (h.ref_abs_in-h.ref_evap_out)/(h.RHEXideal-h.ref_evap_out)
    eta.SHEX = (h.sol_des_in-h.sol_pump_out)/(h.sol_des_out-h.sol_pump_out)
    eta.SHEXI = (h.sol_des_inI-h.sol_pump_outI)/(h.sol_des_outI-h.sol_pump_outI)
    # COP
    PP.COP = Q.evap/(Q.des + PP.W_pump + PP.W_pumpI)
    # "Umlauf"
    PP.f = (m.sol_rich+m.sol_richI)/(m.ref+m.refI)
    # Energy balance
    PP.energyBalance = Q.des + Q.des_mid + Q.evap + PP.W_pump + PP.W_pumpI + Q.cond + Q.condI + Q.abs
    # Mass balance (Absorber)
    PP.massBalance = (m.ref+m.refI) + (m.sol_poor+m.sol_poorI) - (m.sol_rich+m.sol_richI)
    # ----------------------------------------------------------------------- #
    ## Check
    # Refrigerant concentrations
    if (w.NH3_rich < w.NH3_poor):
        print("w_NH3_rich < w_NH3_poor")
        #error("w_NH3_rich < w_NH3_poor")
 
    if (w.NH3_rich < w.NH3_poorI):
        print("w_NH3_rich < w_NH3_poorI")
        #error("w_NH3_rich < w_NH3_poorI")
    
    if (w.NH3_rich < 0):
        print("w_NH3_rich < 0")
        #error("w_NH3_rich < 0")
    
    if (w.NH3_poor < 0):
        print("w_NH3_poor < 0")
        #error("w_NH3_poor < 0")
    
    if (w.NH3_poorI < 0):
        print("w_NH3_poorI < 0")
        #error("w_NH3_poorI < 0")
    
    if (w.NH3_rich - w.NH3_poor < 0.005):
        print("w_NH3_rich - w_NH3_poor < 0.005")
        #    error("w_NH3_rich - w_NH3_poor < 0.005")
    
    if (w.NH3_rich - w.NH3_poorI < 0.005):
        print("w_NH3_rich - w_NH3_strongI < 0.005")
        #    error("w_NH3_rich - w_NH3_strongI < 0.005")
    
    # Mass flows
    if (m.ref<0 or m.sol_poor<0 or m.sol_rich<0):
        print("mass flow is negativ")
        #error("mass flow is negativ")
    
    if (m.refI<0 or m.sol_poorI<0 or m.sol_richI<0):
        print("mass flow I is negativ")
        #error("mass flow I is negativ")
    
    # Energy and mass balance
    if (abs(PP.energyBalance) > 1):
        print("Energy is not conserved")
        #error("Energy is not conserved")
    
    if (abs(PP.massBalance) > 1):
        print("Mass is not conserved")
        #error("Mass is not conserved")
    
    #-------------------------------------------------------------------------#
    
    return T, p, h, m, w, eta, Q, PP, s

