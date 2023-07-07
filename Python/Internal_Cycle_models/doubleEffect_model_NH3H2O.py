import CoolProp.CoolProp as CoolProp
import math
import numpy
from Fluids import LiBrSol, NH3H2O
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import os
import sys

os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])


def doubleEffect_model_NH3H2O(T, p, h, m, eta, Q, HX, s):
    """
    # Function Double Effect Model NH3H2O (Parallel Flow Configuration)
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

    # Input/ Output
    Input:
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
    - HX.dT_ref_des
    - HX.dT_ref_desI

    Output:
    1. Temperature struct                       [K]
    2. Pressure struct                          [Pa]
    3. Specific enthalpy struct                 [J/kg]
    4. Mass flow rate struct                    [kg/s]
    5. Mass fraction struct                     [kg/s]
    6. Efficiency struct                        [-]
    7. Heat flow struct                         [J]
    8. Post process struct
    9. Setup struct

    # Absorption System
    Components:
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

    # Assumptions
    - Temperature of refrigerant leaving desorber is 5K below desorber temp.
    - All components of the system operate in steady state
    - Solution leaving desorber and absorber is saturated
    - Refrigerant leaving condenser and evaporator is saturated
    - Pressure drops in the system components are negelcted
    - Heat capacity of solution assumed to be constant in SHEX
    """
    # ----------------------------------------------------------------------- #
    # # Definition of constants
    # ------------------------Necessary Constants---------------------------- #
    M_NH3 = 0.017031                 # [kg/mol]
    M_H2O = 0.018015268              # [kg/mol]
    # ----------------------------------------------------------------------- #
    # # Initialise required variables
    # ----------------------------------------------------------------------- #
    class var:
        i = 0
    w= var()
    cp=var()
    PP = var()
    # ---------------------------CALCULATION--------------------------------- #
    # # Calculate pressures
    p.evap = CoolProp.PropsSI('P', 'T', T.evap, 'Q', 1, 'Ammonia')
    p.cond = CoolProp.PropsSI('P', 'T', T.cond, 'Q', 0, 'Ammonia')
    p.cond_int = CoolProp.PropsSI('P', 'T', T.cond_int, 'Q', 0, 'Ammonia')
    # ----------------------------------------------------------------------- #
    # # Find minimal possible internal condensation temperature
    w.NH3_rich = NH3H2O.NH3inSolution_Calc_X_PT_REFPROP(p.evap, T.sol_abs_out)
    w.NH3_poor = NH3H2O.NH3inSolution_Calc_X_PT_REFPROP(p.cond, T.sol_des_out)
    if w.NH3_rich < w.NH3_poor + 0.05:
        w.NH3_poor = w.NH3_rich - 0.05 # Minimal solution concentration difference of 5 %
        T.sol_des_out = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="PQ", hOut="T", iUnits=2, iMass=1, iFlag=0, a=p.cond*10**(-6), b=0, z=[w.NH3_poor, (1-w.NH3_poor)]).Output[0]
        T.cond_int = T.sol_des_out + HX.T_PP_cond_int
        p.cond_int = CoolProp.PropsSI('P', 'T', T.cond_int, 'Q', 0, 'Ammonia')
    # ----------------------------------------------------------------------- #
    # # Refrigerant line
    # Desorber
    T.ref_des_outI = T.sol_des_outI - HX.dT_ref_des
    h.ref_des_outI = CoolProp.PropsSI('H', 'T', T.ref_des_outI, 'P', p.cond_int, 'Ammonia')
    T.ref_des_out = T.sol_des_out - HX.dT_ref_desI
    h.ref_des_out = CoolProp.PropsSI('H', 'T', T.ref_des_out, 'P', p.cond, 'Ammonia')
    # Condenser
    T.ref_cond_inI = T.ref_des_outI
    h.ref_cond_inI = h.ref_des_outI
    T.ref_cond_outI = T.cond_int
    h.ref_cond_outI = CoolProp.PropsSI('H', 'P', p.cond_int, 'Q', 0, 'Ammonia')
    T.ref_cond_in = T.ref_des_out
    h.ref_cond_in = h.ref_des_out
    if T.ext_cond_in + HX.T_PP_cond < T.cond:
        T.ref_cond_out = T.ext_cond_in + HX.T_PP_cond  # Sub-cooling as low as possible
        h.ref_cond_out = CoolProp.PropsSI('H', 'T', T.ref_cond_out, 'P', p.cond, 'Ammonia')
    else:
        T.ref_cond_out = T.cond
        h.ref_cond_out = CoolProp.PropsSI('H', 'T', T.ref_cond_out, 'Q', 0, 'Ammonia')
    # Evaporator
    T.ref_evap_out = T.evap
    h.ref_evap_out = CoolProp.PropsSI('H', 'P', p.evap, 'Q', 1, 'Ammonia')
    # Sub-cooler
    if T.ref_cond_out - HX.T_PP_RHEX > T.ref_evap_out:
        T.ref_abs_in = T.ref_cond_out - HX.T_PP_RHEX
        h.ref_abs_in = CoolProp.PropsSI('H', 'T', T.ref_abs_in, 'P', p.evap, 'Ammonia')
    else:
        T.ref_abs_in = T.ref_evap_out
        h.ref_abs_in = h.ref_evap_out
    h.ref_valve_in = h.ref_cond_out - (h.ref_abs_in - h.ref_evap_out)
    T.ref_valve_in = CoolProp.PropsSI('T', 'H', h.ref_valve_in, 'P', p.cond, 'Ammonia')
    # Expansion valve
    h.ref_valve_outI = h.ref_cond_outI  # Isenthalpic thottle
    T.ref_valve_outI = CoolProp.PropsSI('T', 'H', h.ref_valve_outI, 'P', p.cond, 'Ammonia')
    h.ref_evap_in = h.ref_valve_in  # Isenthalpic thottle
    T.ref_evap_in = CoolProp.PropsSI('T', 'H', h.ref_evap_in, 'P', p.evap, 'Ammonia')
    # ----------------------------------------------------------------------- #
    # # Rich solution (High ref. concentration)
    # Absorber
    w.NH3_rich = NH3H2O.NH3inSolution_Calc_X_PT(p.evap, T.sol_abs_out)
    h.sol_abs_out = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="PT", hOut="H", iUnits=2,iMass=1, iFlag=0, a=p.evap*10**(-6), b=T.sol_abs_out, z=[w.NH3_rich, (1-w.NH3_rich)]).Output[0] * 1000
    # Pump
    # Low pressure
    roh_sol_abs_out = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="TQ", hOut="D", iUnits=2, iMass=1, iFlag=0, a=T.sol_abs_out, b=0, z=[w.NH3_rich, (1-w.NH3_rich)]).Output[0]
    h.sol_pump_out = h.sol_abs_out + (1 / roh_sol_abs_out) * (p.cond - p.evap) / eta.pump
    T.sol_pump_out = T.sol_abs_out  # only small temperature change
    cp.sol_pump_out = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="PH", hOut="C", iUnits=2,iMass=1, iFlag=0, a=p.cond*10**(-6), b=h.sol_pump_out, z=[w.NH3_rich, (1-w.NH3_rich)]).Output[0] * 1000
    # High pressure
    h.sol_pump_outI = h.sol_abs_out + (1 / roh_sol_abs_out) * (p.cond_int - p.evap) / eta.pump
    T.sol_pump_outI = T.sol_abs_out  # only small temperature change
    cp.sol_pump_outI = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="PH", hOut="C", iUnits=2,iMass=1, iFlag=0, a=p.cond_int*10**(-6), b=h.sol_pump_outI, z=[w.NH3_rich, (1-w.NH3_rich)]).Output[0] * 1000
    # ----------------------------------------------------------------------- #
    # # Poor solution (Low ref. concentration)
    # Low pressure
    # Desorber
    w.NH3_poor = NH3H2O.NH3inSolution_Calc_X_PT_REFPROP(p.cond, T.sol_des_out)
    h.sol_des_out = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="TQ", hOut="H", iUnits=2, iMass=1, iFlag=0, a=T.sol_des_out, b=0, z=[w.NH3_poor, (1-w.NH3_poor)]).Output[0] * 1000
    # SHEX
    if T.sol_pump_out + HX.T_PP_SHEX < T.sol_des_out:
        T.sol_valve_in = T.sol_pump_out + HX.T_PP_SHEX
    else:
        T.sol_valve_in = T.sol_des_out
    h.sol_valve_in = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="TP", hOut="H", iUnits=2, iMass=1, iFlag=0, a=T.sol_valve_in, b=p.cond*10**(-6), z=[w.NH3_poor, (1-w.NH3_poor)]).Output[0] * 1000
    # High pressure
    # DesorberI
    w.NH3_poorI = NH3H2O.NH3inSolution_Calc_X_PT_REFPROP(p.cond_int, T.sol_des_outI)
    h.sol_des_outI = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="TQ", hOut="H", iUnits=2, iMass=1, iFlag=0, a=T.sol_des_outI, b=0, z=[w.NH3_poorI, (1-w.NH3_poorI)]).Output[0] * 1000
    # SHEXI
    if T.sol_pump_outI + HX.T_PP_SHEXI < T.sol_des_outI:
        T.sol_valve_inI = T.sol_pump_outI + HX.T_PP_SHEXI
    else:
        T.sol_valve_inI = T.sol_des_outI
    h.sol_valve_inI = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="TP", hOut="H", iUnits=2, iMass=1, iFlag=0, a=T.sol_valve_inI, b=p.cond_int*10**(-6), z=[w.NH3_poorI, (1-w.NH3_poorI)]).Output[0] * 1000
    # ----------------------------------------------------------------------- #
    # # Energy balances and mass conservation
    # Solve system of linear equations (energy conservation, mass conservation)
    # Solution vector: Q.condI, Q.des_mid, m.ref, m.sol_rich, m.sol_poor, m.refI, m.sol_richI, m.sol_poorI
    match s.requirement:
        case 'Q_des':
            A = numpy.array([[   0,          0,  0,                  0,              0,              -h.ref_des_out,                     h.sol_pump_outI,    -h.sol_valve_inI],
                                [0,          0,  0,                  0,              0,              -1,                                 w.NH3_rich,         -w.NH3_poorI],
                                [0,          0,  0,                  0,              0,              -1,                                 1,                  -1],
                                [1,          0,  0,                  0,              0,              (h.ref_cond_inI-h.ref_cond_outI),   0,                  0],
                                [1,          1,  0,                  0,              0,              0,                                  0,                  0],
                                [0,          1,  -h.ref_des_out,     h.sol_pump_out, -h.sol_valve_in,0,                                  0,                  0],
                                [0,          0,  -1,                 w.NH3_rich,     -w.NH3_poor,    0,                                  0,                  0],
                                [0,          0,  -1,                 1,              -1,             0,                                  0,                  0]])
            b = numpy.array([    -Q.dec,     0,  0,                  0,              0,              0,                                  0,                  0])
        case 'Q_evap':
            A = numpy.array([  [ 0,          0,  (h.ref_evap_out-h.ref_evap_in),     0,              0,                  (h.ref_evap_out-h.ref_evap_in),     0,              0],
                                [0,          0,  0,                                  0,              0,                  -1,                                 w.NH3_rich,    -w.NH3_poorI],
                                [0,          0,  0,                                  0,              0,                  -1,                                 1,              -1],
                                [1,          0,  0,                                  0,              0,                  (h.ref_cond_inI-h.ref_cond_outI),   0,              0],
                                [1,          1,  0,                                  0,              0,                  0,                                  0,              0],
                                [0,          1,  -h.ref_des_out,                     h.sol_pump_out, -h.sol_valve_in,    0,                                  0,              0],
                                [0,          0,  -1,                                 w.NH3_rich,     -w.NH3_poor,        0,                                  0,              0],
                                [0,          0,  -1,                                 1,              -1,                 0,                                  0,              0]])
            b = numpy.array([    Q.dec,      0,  0,                                  0,              0,                  0,                                  0,              0])
        case _:
            sys.exit('AKM decider is not defined properly. Use Q_des or Q_evap')
    y = numpy.linalg.solve(A, b)
    # High pressure side
    Q.condI = y[0]
    m.refI = y[5]
    m.sol_richI = y[6]
    m.sol_poorI = y[7]
    # Low pressure side
    Q.des_mid = y[1]
    m.ref = y[2]
    m.sol_rich = y[3]
    m.sol_poor = y[4]
    # ----------------------------------------------------------------------- #
    # # Rich solution after SHEX
    # Low pressure
    h.sol_des_in = (m.sol_poor*h.sol_des_out + m.sol_rich*h.sol_pump_out - m.sol_poor*h.sol_valve_in) / m.sol_rich
    T.sol_des_in = T.sol_pump_out + (h.sol_des_in - h.sol_pump_out)/cp.sol_pump_out
    # High pressure
    h.sol_des_inI = (m.sol_poorI*h.sol_des_outI + m.sol_richI*h.sol_pump_outI - m.sol_poorI*h.sol_valve_inI) / m.sol_richI
    T.sol_des_inI = T.sol_pump_outI + (h.sol_des_inI - h.sol_pump_outI)/cp.sol_pump_outI
    # # Poor solution after valve
    # Low pressure
    h.sol_abs_in = h.sol_valve_in
    T.sol_abs_in = T.sol_valve_in
    # High pressure
    h.sol_abs_inI = h.sol_valve_inI
    T.sol_abs_inI = T.sol_valve_inI
    # ----------------------------------------------------------------------- #
    # # Check
    # Solution concentrations
    if w.NH3_rich < 0:
        sys.exit("w_NH3_rich < 0")
    if w.NH3_poor < 0:
        sys.exit("w_NH3_poor < 0")
    if w.NH3_rich - w.NH3_poor < 0.005:
        sys.exit("w_NH3_rich - w_NH3_poor < 0.005")
    if w.NH3_poorI < 0:
        sys.exit("w_NH3_poorI < 0")
    if w.NH3_rich < w.NH3_poorI:
        sys.exit("w_NH3_rich < w_NH3_poorI")
    if w.NH3_rich - w.NH3_poorI < 0.005:
        sys.exit("w_NH3_rich - w_NH3_strongI < 0.005")
    # Mass flows
    if m.ref < 0 or m.sol_poor < 0 or m.sol_rich < 0:
        sys.exit("mass flow is negativ")
    if m.refI < 0 or m.sol_poorI < 0 or m.sol_richI < 0:
        sys.exit("mass flow I is negativ")
    # ----------------------------------------------------------------------- #
    # # Post processing
    # Calculate fluxes over system boundary
    Q.cond = h.ref_cond_out * (m.ref + m.refI) - h.ref_cond_in * m.ref - h.ref_valve_outI * m.refI
    Q.evap = (m.ref + m.refI) * (h.ref_evap_out - h.ref_evap_in)
    Q.abs = (m.sol_rich + m.sol_richI) * h.sol_abs_out - (m.ref + m.refI) * h.ref_abs_in - m.sol_poor * h.sol_abs_in - m.sol_poorI * h.sol_abs_inI
    Q.des = m.refI * h.ref_des_outI + m.sol_poorI * h.sol_valve_inI - m.sol_richI * h.sol_pump_outI
    PP.W_pump = m.sol_rich * (h.sol_pump_out - h.sol_abs_out)
    PP.W_pumpI = m.sol_richI * (h.sol_pump_outI - h.sol_abs_out)
    # Heat Exchanger
    Q.SHEX = m.sol_poor * (h.sol_des_out - h.sol_valve_in)
    Q.SHEXI = m.sol_poorI * (h.sol_des_outI - h.sol_valve_inI)
    Q.RHEX = m.ref * (h.ref_abs_in - h.ref_evap_out)
    h.RHEXideal = CoolProp.PropsSI('H', 'T', T.ref_cond_out, 'P', p.evap, 'Ammonia')
    eta.RHEX = (h.ref_abs_in - h.ref_evap_out) / (h.RHEXideal - h.ref_evap_out)
    eta.SHEX = (h.sol_des_in - h.sol_pump_out) / (h.sol_des_out - h.sol_pump_out)
    eta.SHEXI = (h.sol_des_inI - h.sol_pump_outI) / (h.sol_des_outI - h.sol_pump_outI)
    # COP
    PP.COP = Q.evap / (Q.des + PP.W_pump + PP.W_pumpI)
    # "Umlauf"
    PP.f = (m.sol_rich + m.sol_richI) / (m.ref + m.refI)
    # Energy balance
    PP.energyBalance = Q.des + Q.des_mid + Q.evap + PP.W_pump + PP.W_pumpI + Q.cond + Q.condI + Q.abs
    # Mass balance (Absorber)
    PP.massBalance = (m.ref + m.refI) + (m.sol_poor + m.sol_poorI) - (m.sol_rich + m.sol_richI)
    # ----------------------------------------------------------------------- #    
    # Energy and mass balance
    if abs(PP.energyBalance) > 1:
        sys.exit("Energy is not conserved")
    if abs(PP.massBalance) > 1:
        sys.exit("Mass is not conserved")
    # ----------------------------------------------------------------------- #
    return T, p, h, m, w, eta, Q, PP, s
