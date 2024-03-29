import CoolProp.CoolProp as CoolProp
from Fluids import NH3H2O
import numpy
import math
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import os
import sys

# Initialize REFPROP
os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])


def doubleLift_model_NH3H2O(T, p, h, m, eta, Q, HX, s):
    """
    # Function Double Lift Model NH3H2O
    Author  : Ludwig Irrgang
    Date    : 28.07.2022
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
    - 2 x Absorber
    - 2 x Pump
    - 3 x Throttle Valves
    - Condenser
    - Evaporator
    - Solution Heat Exchanger
    - Refrigerant Heat Exchanger
    - Working fluid: Ammonia Water Solution
    - Refrigerant: Ammonia

    # Assumptions
    - Temperature of refrigerant leaving desorber is 5K below desorber temp.
    - All components of the system operate in steady state
    - Solution leaving generator and absorber is saturated
    - Refrigerant leaving condenser and evaporator is saturated
    - Pressure drops in the system components are neglected
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
    w = var()
    cp = var()
    PP = var()
    # ---------------------------CALCULATION--------------------------------- #
    # # Calculate pressures
    p.evap = CoolProp.PropsSI('P', 'T', T.evap, 'Q', 1, 'AMMONIA')
    p.cond = CoolProp.PropsSI('P', 'T', T.cond, 'Q', 0, 'AMMONIA')
    p.mid = math.sqrt(p.evap*p.cond)  # Venegas
    # ----------------------------------------------------------------------- #
    # # Refrigerant line
    # Desorber
    T.ref_des_out = T.sol_des_out - HX.dT_ref_des
    h.ref_des_out = CoolProp.PropsSI('H', 'T', T.ref_des_out, 'P', p.mid, 'AMMONIA')
    T.ref_des_outI = T.sol_des_outI - HX.dT_ref_desI
    h.ref_des_outI = CoolProp.PropsSI('H', 'T', T.ref_des_outI, 'P', p.cond, 'AMMONIA')
    # Condenser
    T.ref_cond_in = T.ref_des_outI
    h.ref_cond_in = h.ref_des_outI
    if T.ext_cond_in + HX.T_PP_cond < T.cond:
        T.ref_cond_out = T.ext_cond_in + HX.T_PP_cond  # Sub-cooling as low as possible
        h.ref_cond_out = CoolProp.PropsSI('H', 'T', T.ref_cond_out, 'P', p.cond, 'AMMONIA')
    else:
        T.ref_cond_out = T.cond
        h.ref_cond_out = CoolProp.PropsSI('H', 'T', T.ref_cond_out, 'Q', 0, 'AMMONIA')
    # Evaporator
    T.ref_evap_out = T.evap
    h.ref_evap_out = CoolProp.PropsSI('H', 'P', p.evap, 'Q', 1, 'AMMONIA')
    # Sub-cooler (heat capacity of steam lower than liquid)
    if T.ref_cond_out - HX.T_PP_RHEX > T.ref_evap_out:
        T.ref_abs_in = T.ref_cond_out - HX.T_PP_RHEX
        h.ref_abs_in = CoolProp.PropsSI('H', 'T', T.ref_abs_in, 'P', p.evap, 'AMMONIA')
    else:
        T.ref_abs_in = T.ref_evap_out
        h.ref_abs_in = h.ref_evap_out
    h.ref_valve_in = h.ref_cond_out - (h.ref_abs_in - h.ref_evap_out)
    T.ref_valve_in = CoolProp.PropsSI('T', 'H', h.ref_valve_in, 'P', p.cond, 'AMMONIA')
    # Throttle
    h.ref_evap_in = h.ref_valve_in  # Isenthalpic throttle
    T.ref_evap_in = CoolProp.PropsSI('T', 'H', h.ref_evap_in, 'P', p.evap, 'AMMONIA')
    # ----------------------------------------------------------------------- #
    # # Rich solution (High ref. concentration)
    # Low pressure cycle
    # Absorber
    w.NH3_rich = NH3H2O.NH3inSolution_Calc_X_PT_REFPROP(p.evap, T.sol_abs_out)
    h.sol_abs_out = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="TP", hOut="H", iUnits=2, iMass=1, iFlag=0, a=T.sol_abs_out, b=p.evap*10**(-6), z=[w.NH3_rich, (1-w.NH3_rich)]).Output[0] * 1000
    # Pump
    roh_sol_abs_out = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="TQ", hOut="D", iUnits=2, iMass=1, iFlag=0, a=T.sol_abs_out, b=0, z=[w.NH3_rich, (1-w.NH3_rich)]).Output[0]
    h.sol_pump_out = h.sol_abs_out + (1 / roh_sol_abs_out) * (p.mid - p.evap) / eta.pump
    T.sol_pump_out = T.sol_abs_out  # only small temperature change
    cp.sol_pump_out = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="PH", hOut="C", iUnits=2, iMass=1, iFlag=0, a=p.mid*10**(-6), b=h.sol_pump_out, z=[w.NH3_rich, (1-w.NH3_rich)]).Output[0] * 1000
    # High pressure cycle
    # Absorber
    w.NH3_richI = NH3H2O.NH3inSolution_Calc_X_PT_REFPROP(p.mid,T.sol_abs_outI)
    h.sol_abs_outI = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="TP", hOut="H", iUnits=2, iMass=1, iFlag=0, a=T.sol_abs_outI, b=p.mid*10**(-6), z=[w.NH3_richI, (1-w.NH3_richI)]).Output[0] * 1000
    # Pump
    roh_sol_abs_outI = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="TQ", hOut="D", iUnits=2, iMass=1, iFlag=0, a=T.sol_abs_outI, b=0, z=[w.NH3_richI, (1-w.NH3_richI)]).Output[0]
    h.sol_pump_outI = h.sol_abs_outI + (1 / roh_sol_abs_outI) * (p.cond - p.mid) / eta.pump
    T.sol_pump_outI = T.sol_abs_outI  # only small temperature change
    cp.sol_pump_outI = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="PH", hOut="C", iUnits=2, iMass=1, iFlag=0, a=p.mid*10**(-6), b=h.sol_pump_outI, z=[w.NH3_richI, (1-w.NH3_richI)]).Output[0] * 1000
    # ----------------------------------------------------------------------- #
    # # Poor solution (Low ref. concentration)
    # High pressure
    # Desorber
    w.NH3_poorI = NH3H2O.NH3inSolution_Calc_X_PT_REFPROP(p.cond, T.sol_des_outI)
    h.sol_des_outI = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="TQ", hOut="H", iUnits=2, iMass=1, iFlag=0, a=T.sol_des_outI, b=0, z=[w.NH3_poorI, (1 - w.NH3_poorI)]).Output[0] * 1000
    # SHEX
    if T.sol_pump_outI + HX.T_PP_SHEXI < T.sol_des_outI:
        T.sol_valve_inI = T.sol_pump_outI + HX.T_PP_SHEXI
    else:
        T.sol_valve_inI = T.sol_des_outI
    h.sol_valve_inI = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="TP", hOut= "H", iUnits=2, iMass=1, iFlag=0, a=T.sol_valve_inI, b=p.cond*10**(-6), z=[w.NH3_poorI, (1-w.NH3_poorI)]).Output[0] * 1000
    # Low pressure
    # Poor solution leaving desorber
    w.NH3_poor = NH3H2O.NH3inSolution_Calc_X_PT_REFPROP(p.mid, T.sol_des_out)
    h.sol_des_out = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="TQ", hOut="H", iUnits=2, iMass=1, iFlag=0, a=T.sol_des_out, b=0, z=[w.NH3_poor, (1-w.NH3_poor)]).Output[0] * 1000
    # SHEX
    if T.sol_pump_out + HX.T_PP_SHEX < T.sol_des_out:
        T.sol_valve_in = T.sol_pump_out + HX.T_PP_SHEX
    else:
        T.sol_valve_in = T.sol_des_out
    h.sol_valve_in = RP.REFPROPdll(hFld="AMMONIA;WATER", hIn="TP", hOut="H", iUnits=2, iMass=1, iFlag=0, a=T.sol_valve_in, b=p.mid*10**(-6), z=[w.NH3_poor, (1-w.NH3_poor)]).Output[0] * 1000
    # ----------------------------------------------------------------------- #
    # # Energy balances and mass conservation
    # Solve system of linear equations (energy conservation, mass conservation)
    # Solution vector: Q.des, m.ref, m.sol_rich, m.sol_poor, Q.desI, m.refI, m.sol_richI, m.sol_poorI
    match s.requirement:
        case 'Q_des':
            A = numpy.array([[   1,      0,              0,              0,                  1,  0,                  0,                  0],
                                [0,      -1,             1,              -1,                 0,  0,                  0,                  0],
                                [1,      -h.ref_des_out, h.sol_pump_out, -h.sol_valve_in,    0,  0,                  0,                  0],
                                [0,      -1,             w.NH3_rich,     -w.NH3_poor,        0,  0,                  0,                  0],
                                [0,      0,              0,              0,                  0,  -1,                 1,                  -1],
                                [0,      0,              0,              0,                  1,  -h.ref_des_outI,    h.sol_pump_outI,    -h.sol_valve_inI],
                                [0,      0,              0,              0,                  0,  -1,                 w.NH3_richI,        -w.NH3_poorI],
                                [0,      1,              0,              0,                  0,  -1,                 0,                  0]])
            b = numpy.array([    Q.dec,  0,              0,              0,                  0,  0,                  0,                  0])
        case 'Q_evap':
            A = numpy.array([[   1,      0,              0,              0,                  0,  (h.ref_evap_out-h.ref_evap_in), 0,                  0],
                                [0,      -1,             1,              -1,                 0,  0,                              0,                  0],
                                [1,      -h.ref_des_out, h.sol_pump_out, -h.sol_valve_in,    0,  0,                              0,                  0],
                                [0,      -1,             w.NH3_rich,     -w.NH3_poor,        0,  0,                              0,                  0],
                                [0,      0,              0,              0,                  0,  -1,                             1,                  -1],
                                [0,      0,              0,              0,                  1,  -h.ref_des_outI,                h.sol_pump_outI,    -h.sol_valve_inI],
                                [0,      0,              0,              0,                  0,  -1,                             w.NH3_richI,        -w.NH3_poorI],
                                [0,      1,              0,              0,                  0,  -1,                             0,                  0]])
            b = numpy.array([    Q.dec,  0,              0,              0,                  0,  0,                              0,                  0])
        case _:
            sys.exit('AKM decider is not defined properly. Use Q_des or Q_evap')
    y = numpy.linalg.solve(A, b)
    # Low pressure side
    Q.des = y[0]
    m.ref = y[1]
    m.sol_rich = y[2]
    m.sol_poor = y[3]
    # High pressure side
    Q.desI = y[4]
    m.refI = y[5]
    m.sol_richI = y[6]
    m.sol_poorI = y[7]
    # ----------------------------------------------------------------------- #
    # # Rich solution after SHEX
    # Low pressure
    h.sol_des_in = (m.sol_poor * h.sol_des_out + m.sol_rich * h.sol_pump_out - m.sol_poor * h.sol_valve_in) / m.sol_rich
    T.sol_des_in = T.sol_pump_out + (h.sol_des_in - h.sol_pump_out) / cp.sol_pump_out
    # High pressure
    h.sol_des_inI = (m.sol_poorI * h.sol_des_outI + m.sol_richI * h.sol_pump_outI - m.sol_poorI * h.sol_valve_inI) / m.sol_richI
    T.sol_des_inI = T.sol_pump_outI + (h.sol_des_inI - h.sol_pump_outI) / cp.sol_pump_outI
    # # Poor solution after valve
    # Low pressure
    h.sol_abs_in = h.sol_valve_in
    T.sol_abs_in = T.sol_valve_in
    # High pressure
    h.sol_abs_inI = h.sol_valve_inI
    T.sol_abs_inI = T.sol_valve_inI
    # ----------------------------------------------------------------------- #
    # # Check
    if w.NH3_rich < 0:
        sys.exit("w_NH3_rich < 0")
    if w.NH3_richI < 0:
        sys.exit("w_NH3_richI < 0")
    if w.NH3_richI - w.NH3_poorI < 0.005:
        sys.exit("w_NH3_richI - w_NH3_poorI < 0.005")
    if w.NH3_poorI < 0:
        sys.exit("w_NH3_poorI < 0")
    if w.NH3_rich < w.NH3_poor:
        sys.exit("w_NH3_rich < w_NH3_poor")
    if w.NH3_poor < 0:
        sys.exit("w_NH3_poor < 0")
    if w.NH3_rich - w.NH3_poor < 0.005:
        sys.exit("w_NH3_rich - w_NH3_poor < 0.001")
    # Mass flows
    if m.ref < 0 or m.sol_poor < 0 or m.sol_rich < 0:
        sys.exit("mass flow is negativ")
    if m.refI < 0 or m.sol_poorI < 0 or m.sol_richI < 0:
        sys.exit("mass flow I is negativ")
    # ----------------------------------------------------------------------- #
    # # Post processing
    # Fluxes over system boundary
    Q.cond = m.refI * (h.ref_cond_out - h.ref_cond_in)
    Q.evap = m.refI * (h.ref_evap_out - h.ref_evap_in)
    Q.abs = m.sol_rich * h.sol_abs_out - m.refI * h.ref_abs_in - m.sol_poor * h.sol_abs_in
    Q.absI = m.sol_richI * h.sol_abs_outI - m.ref * h.ref_des_out - m.sol_poorI * h.sol_abs_inI
    PP.W_pump = m.sol_rich * (h.sol_pump_out - h.sol_abs_out)
    PP.W_pumpI = m.sol_richI * (h.sol_pump_outI - h.sol_abs_outI)
    Q.des = m.ref * h.ref_des_out + m.sol_poor * h.sol_des_out - m.sol_rich * h.sol_des_in
    Q.desI = m.refI * h.ref_des_outI + m.sol_poorI * h.sol_des_outI - m.sol_richI * h.sol_des_inI
    Q.des_ges = Q.des + Q.desI
    # Heat Exchanger
    Q.SHEX = m.sol_poor * (h.sol_des_out - h.sol_valve_in)
    Q.SHEXI = m.sol_poorI * (h.sol_des_outI - h.sol_valve_inI)
    Q.RHEX = m.refI * (h.ref_abs_in - h.ref_evap_out)
    h.RHEXideal = CoolProp.PropsSI('H', 'T', T.ref_cond_out, 'P', p.evap, 'AMMONIA')
    eta.RHEX = (h.ref_abs_in - h.ref_evap_out) / (h.RHEXideal - h.ref_evap_out)
    eta.SHEX = (h.sol_des_in - h.sol_pump_out) / (h.sol_des_out - h.sol_pump_out)
    eta.SHEXI = (h.sol_des_inI - h.sol_pump_outI) / (h.sol_des_outI - h.sol_pump_outI)
    # COP
    PP.COP = Q.evap / (Q.des + Q.desI + PP.W_pump + PP.W_pumpI)
    # "Umlauf"
    PP.f = (m.sol_rich + m.sol_richI) / m.ref
    # Energy balance
    PP.energyBalance = Q.des + Q.desI + Q.evap + PP.W_pump + PP.W_pumpI + Q.cond + Q.abs + Q.absI
    # Mass balance
    PP.massBalance = m.ref + m.sol_poor - m.sol_rich + m.refI + m.sol_poorI - m.sol_richI
    # ----------------------------------------------------------------------- #
    # # Check
    # Energy and mass balance
    if abs(PP.energyBalance) > 1:
        sys.exit("Energy is not conserved")
    if abs(PP.massBalance) > 1:
        sys.exit("Mass is not conserved")
    # ----------------------------------------------------------------------- #
    return T, p, h, m,w, eta, Q, PP, s
