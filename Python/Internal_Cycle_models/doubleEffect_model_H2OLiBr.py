import CoolProp.CoolProp as CoolProp
import math
import numpy
from Fluids import LiBrSol, NH3H2O, H2O


def doubleEffect_model_H2OLiBr(T, p, h, m, eta, Q, HX, s):
    """
    ## Function Double Effect Model H2OLiBr (Parallel Flow Configuration)
    # ----------------------------------------------------------------------- #
    #{
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
    #}
    # ----------------------------------------------------------------------- #
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
    1. Temperature Struct                       [K]
    2. Pressure Struct                          [Pa]
    3. Specific enthalpy Struct                 [J/kg]
    4. Mass flow rate Struct                    [kg/s]
    5. Mass fraction Struct                     [kg/s]
    6. Efficiency Struct                        [-]
    7. Heat flow Struct                         [J]
    8. Post Process Struct
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
    - 3 x Throttle valves
    - Evaporator
    - Absorber
    - Solution Heat exchanger
    - Refrigerant Heat exchanger
    - Working fluid: Water lithium-bromide solution
    - Refrigerant: Water
    #}
    # ----------------------------------------------------------------------- #
    ## Annotations
    # ----------------------------------------------------------------------- #
    #{ 
    - x stands for molar fraction
    - w stands for mass fraction
    - Functions for LiBr are using molar state properties and molar fraction
    #} 
    # ----------------------------------------------------------------------- #
    ## Assumptions
    # ----------------------------------------------------------------------- #
    #{ 
    - Temperature of refrigerant leaving desorber is 5K below desorber temp.
    - All components of the system operate in steady state
    - Solution leaving generator and absorber is saturated
    - Refrigerant leaving evaporator is saturated
    - Pressure drops in the system components are negelcted
    - Heat capacity of solution assumed to be constant in undersaturated solutions
    #} 
    """
    # ----------------------------------------------------------------------- #
    ## Define Variables
    # ------------------------Necessary Constants---------------------------- #
    M_LiBr = 0.08685               #[kg/mol]
    M_H2O = 0.018015268            #[kg/mol]

    # -----------------------------------------------------------------------  #
    # # Initialisation of variable classes for storing calculatioin results
    # -----------------------------------------------------------------------  #
    class var:
        i = 0 
    
    x = var() 
    w = var() 
    rho = var() 
    v = var() 
    cp = var() 
    PP = var()
    ## Calculation
    # ---------------------------CALCULATION--------------------------------- #
    ## Calculate pressures
    p.evap = CoolProp.PropsSI('P','T',T.evap,'Q',1,'Water')
    p.cond = CoolProp.PropsSI('P','T',T.cond,'Q',0,'Water')
    p.cond_int = CoolProp.PropsSI('P','T',T.cond_int,'Q',0,'Water')
    # ----------------------------------------------------------------------- #
    ## Rich solution (High ref. concentration)
    # Absorber
    x.LiBr_rich = LiBrSol.Calc_X_from_T_p_satLiBrSol_Patek(T.sol_abs_out,p.evap)
    x.H2O_rich = 1 - x.LiBr_rich
    w.LiBr_rich = x.LiBr_rich*M_LiBr/ (x.LiBr_rich*M_LiBr + x.H2O_rich*M_H2O)
    w.H2O_rich = 1 - w.LiBr_rich
    h.sol_abs_out_mol = LiBrSol.Calc_h_from_T_X_LiBrSol_Patek(T.sol_abs_out,x.LiBr_rich)
    h.sol_abs_out = h.sol_abs_out_mol/(x.LiBr_rich*M_LiBr + x.H2O_rich*M_H2O)
    rho.sol_abs_out_mol = LiBrSol.Calc_rho_from_T_X_LiBrSol_Patek(T.sol_abs_out,x.LiBr_rich)
    rho.sol_abs_out = rho.sol_abs_out_mol*(x.LiBr_rich*M_LiBr + x.H2O_rich*M_H2O)
    v.sol_abs_out = 1/rho.sol_abs_out
    # ----------------------------------------------------------------------- #
    # Prohibit crystallization in SHEXI - necessary internal condensation pressure increase 
    T.sol_pump_outI = T.sol_abs_out # Is assumed to be isothermal
    # SHEX
    if(T.sol_pump_outI + HX.T_PP_SHEXI < T.sol_des_outI):
        T.sol_valve_inI = T.sol_pump_outI + HX.T_PP_SHEXI
    else:
        T.sol_valve_inI = T.sol_des_outI

    w.crI = LiBrSol.crystallization_H2OLiBr("T",T.sol_valve_inI)
    x.crI = w.crI/M_LiBr/ (w.crI/M_LiBr+(1-w.crI)/M_H2O)
    p.crI = LiBrSol.Calc_p_from_T_X_LiBrSol_Patek(T.sol_des_outI,x.crI)
    if (p.cond_int<p.crI):
        p.cond_int = p.crI + 10 # Prevent calculation error
        T.cond_int = CoolProp.PropsSI('T','Q',0,'P',p.cond_int,'Water')

    T.sol_des_out = T.cond_int - HX.T_PP_cond_int
    # ----------------------------------------------------------------------- #
    # Prohibit crystallization in SHEX - necessary condensation pressure increase 
    T.sol_pump_out = T.sol_abs_out # Is assumed to be isothermal
    # SHEX
    if(T.sol_pump_out + HX.T_PP_SHEX < T.sol_des_out):
        T.sol_valve_in = T.sol_pump_out + HX.T_PP_SHEX
    else:
        T.sol_valve_in = T.sol_des_out

    w.cr = LiBrSol.crystallization_H2OLiBr("T",T.sol_valve_in)
    x.cr = w.cr/M_LiBr/ (w.cr/M_LiBr+(1-w.cr)/M_H2O)
    p.cr = LiBrSol.Calc_p_from_T_X_LiBrSol_Patek(T.sol_des_out,x.cr)
    if (p.cond<p.cr):
        p.cond = p.cr + 10 # Prevent calculation error
        T.cond = CoolProp.PropsSI('T','Q',0,'P',p.cond,'Water')
    # ----------------------------------------------------------------------- #
    # Pump
    w.pump = v.sol_abs_out*(p.cond - p.evap)
    h.sol_pump_out = h.sol_abs_out + w.pump/eta.pump
    w.pumpI = v.sol_abs_out*(p.cond_int - p.evap)
    h.sol_pump_outI = h.sol_abs_out + w.pumpI/eta.pump
    #-------------------------------------------------------------------------#
    ## Refrigerant line
    # Desorber
    T.ref_des_out = T.sol_des_out - 5
    h.ref_des_out = CoolProp.PropsSI('H','T',T.ref_des_out,'P',p.cond,'Water')
    T.ref_des_outI = T.sol_des_outI - 5
    h.ref_des_outI = CoolProp.PropsSI('H','T',T.ref_des_outI,'P',p.cond_int,'Water')
    # Condenser
    T.ref_cond_in = T.ref_des_out
    h.ref_cond_in = h.ref_des_out
    T.ref_cond_out = T.ext_cond_in + HX.T_PP_cond # Subcooling as low as possible
    h.ref_cond_out = CoolProp.PropsSI('H','P',p.cond,'T',T.ref_cond_out,'Water')
    T.ref_cond_inI = T.ref_des_outI
    h.ref_cond_inI = h.ref_des_outI
    T.ref_cond_outI = T.cond_int
    h.ref_cond_outI = CoolProp.PropsSI('H','P',p.cond_int,'Q',0,'Water') # Saturated liquid
    # Evaporator
    T.ref_evap_out = T.evap
    h.ref_evap_out = CoolProp.PropsSI('H','P',p.evap,'Q',1,'Water')
    # Subcooler (heat capacity of steam lower than liquid)
    if(T.ref_cond_out - HX.T_PP_RHEX > T.ref_evap_out):
        T.ref_abs_in = T.ref_cond_out - HX.T_PP_RHEX
    else:
        T.ref_abs_in = T.ref_evap_out

    h.ref_abs_in = CoolProp.PropsSI('H','T',T.ref_abs_in,'P',p.evap,'Water') 
    h.ref_valve_in = h.ref_cond_out - (h.ref_abs_in-h.ref_evap_out)
    T.ref_valve_in = CoolProp.PropsSI('T','H',h.ref_valve_in,'P',p.cond,'Water')
    # Expansion valve
    h.ref_valve_outI = h.ref_cond_outI # Isenthalpic thottle
    T.ref_valve_outI = CoolProp.PropsSI('T','H',h.ref_valve_outI,'P',p.cond,'Water')
    h.ref_evap_in = h.ref_valve_in # Isenthalpic thottle
    T.ref_evap_in = CoolProp.PropsSI('T','H',h.ref_evap_in,'P',p.evap,'Water')
    # ----------------------------------------------------------------------- #
    ## Poor solution (Low ref. concentration)
    # Low pressure
    # Desorber
    x.LiBr_poor = LiBrSol.Calc_X_from_T_p_satLiBrSol_Patek(T.sol_des_out,p.cond)
    x.H2O_poor = 1 - x.LiBr_poor
    w.LiBr_poor = x.LiBr_poor*M_LiBr/ (x.LiBr_poor*M_LiBr + x.H2O_poor*M_H2O)
    w.H2O_poor = 1 - w.LiBr_poor
    h.sol_des_out_mol = LiBrSol.Calc_h_from_T_X_LiBrSol_Patek(T.sol_des_out,x.LiBr_poor)
    h.sol_des_out = h.sol_des_out_mol/(x.LiBr_poor*M_LiBr + x.H2O_poor*M_H2O)
    # Assuming constant heat capacity
    cp.sol_des_out_mol = LiBrSol.Calc_cp_from_T_X_LiBrSol_Patek(T.sol_des_out,x.LiBr_poor)
    cp.sol_des_out = cp.sol_des_out_mol/(x.LiBr_poor*M_LiBr + x.H2O_poor*M_H2O)
    h.sol_valve_in = h.sol_des_out+cp.sol_des_out*(T.sol_valve_in-T.sol_des_out)
    # High pressure
    # Desorber
    x.LiBr_poorI = LiBrSol.Calc_X_from_T_p_satLiBrSol_Patek(T.sol_des_outI,p.cond_int)
    x.H2O_poorI = 1 - x.LiBr_poorI
    w.LiBr_poorI = x.LiBr_poorI*M_LiBr/ (x.LiBr_poorI*M_LiBr + x.H2O_poorI*M_H2O)
    w.H2O_poorI = 1 - w.LiBr_poorI
    h.sol_des_out_molI = LiBrSol.Calc_h_from_T_X_LiBrSol_Patek(T.sol_des_outI,x.LiBr_poorI)
    h.sol_des_outI = h.sol_des_out_molI/(x.LiBr_poorI*M_LiBr + x.H2O_poorI*M_H2O)
    # Assuming constant heat capacity
    cp.sol_des_out_molI = LiBrSol.Calc_cp_from_T_X_LiBrSol_Patek(T.sol_des_outI,x.LiBr_poorI)
    cp.sol_des_outI = cp.sol_des_out_molI/(x.LiBr_poorI*M_LiBr + x.H2O_poorI*M_H2O)
    h.sol_valve_inI = h.sol_des_outI+cp.sol_des_outI*(T.sol_valve_inI-T.sol_des_outI)
    # ----------------------------------------------------------------------- #
    ## Energy balances and mass conservation
    # Solve system of linear equations (energy conservation, mass conservation)
    # Solution vector: Q_condI, Q_des_mid, m.ref, m.sol_rich, m.sol_poor, m.refI, m.sol_richI, m.sol_poorI
    match s.requirement:
        case 'Q_des':
            A = numpy.array([[   0,          0,  0,                  0,              0,              -h.ref_des_out,                     h.sol_pump_outI,    -h.sol_valve_inI],
                [    0,          0,  0,                  0,              0,              -1,                                 w.H2O_rich,         -w.H2O_poorI],
                [    0,          0,  0,                  0,              0,              -1,                                 1,                  -1 ],
                [    1,          0,  0,                  0,              0,              (h.ref_cond_inI-h.ref_cond_outI),   0,                  0 ],
                [    1,          1,  0,                  0,              0,              0,                                  0,                  0 ],
                [    0,          1,  -h.ref_des_out,     h.sol_pump_out, -h.sol_valve_in,0,                                  0,                  0 ],
                [    0,          0,  -1,                 w.H2O_rich,     -w.H2O_poor,    0,                                  0,                  0 ],
                [    0,          0,  -1,                 1,              -1,             0,                                  0,                  0   ] ])
            b = numpy.array([   -Q.dec ,    0,  0       ,           0      ,        0,        0,                0,                  0   ])
        case 'Q_evap':
            A = numpy.array([[   0,          0,  (h.ref_evap_out-h.ref_evap_in),     0,              0,                  (h.ref_evap_out-h.ref_evap_in),     0,              0 ],
                [    0,          0,  0,                                  0,              0,                  -1,                                 w.H2O_rich,    -w.H2O_poorI],
                [    0,          0,  0,                                  0,              0,                  -1,                                 1,              -1],
                [    1,          0,  0,                                  0,              0,                  (h.ref_cond_inI-h.ref_cond_outI),   0,              0],
                [    1,          1,  0,                                  0,              0,                  0,                                  0,              0],
                [    0,          1,  -h.ref_des_out,                     h.sol_pump_out, -h.sol_valve_in,    0,                                  0,              0],
                [    0,          0,  -1,                                 w.H2O_rich,     -w.H2O_poor,        0,                                  0,              0],
                [    0,          0,  -1,                                 1,              -1,                 0,                                  0,              0 ]  ])
            b = numpy.array([   Q.dec,      0,  0,                                  0,              0,                  0,                                  0,              0   ])
        case _:
            print('AKM decider is not defined properly. Use Q_des or Q_evap')
            #error('AKM decider is not defined properly. Use Q_des or Q_evap')

    y = numpy.linalg.solve(A,b) # replacement for maltab function "\", mldivide
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
    ## Rich solution after SHEX
    Q.SHEX = (h.sol_des_out - h.sol_valve_in)*m.sol_poor
    h.sol_des_in = h.sol_pump_out + Q.SHEX/m.sol_rich
    T.sat_rich_SHEX = LiBrSol.Calc_T_from_p_X_satLiBrSol_Patek(p.cond,x.LiBr_rich)
    h.sat_rich_SHEX = LiBrSol.Calc_h_from_T_X_LiBrSol_Patek(T.sat_rich_SHEX,x.LiBr_rich)
    if (h.sat_rich_SHEX>h.sol_des_in): # Saturation is not reached in SHEX - no evaporation
        h.sol_des_in = h.sol_des_in
        cp.sol_pump_out_mol = LiBrSol.Calc_cp_from_T_X_LiBrSol_Patek(T.sol_pump_out,x.LiBr_rich)
        cp.sol_pump_out = cp.sol_pump_out_mol/(x.LiBr_rich*M_LiBr + x.H2O_rich*M_H2O)
        T.sol_des_in = T.sol_pump_out + (h.sol_des_in-h.sol_pump_out)/cp.sol_pump_out
    else:
        T.sol_des_in = LiBrSol.Calc_state_SHEX_exit(T.sat_rich_SHEX, p.cond, h.sol_des_in, m.sol_rich, w.H2O_rich, 0.01)

    ## Rich solution after SHEXI
    Q.SHEXI = (h.sol_des_outI - h.sol_valve_inI)*m.sol_poorI
    h.sol_des_inI = h.sol_pump_outI + Q.SHEXI/m.sol_richI
    T.sat_rich_SHEXI = LiBrSol.Calc_T_from_p_X_satLiBrSol_Patek(p.cond_int,x.LiBr_rich)
    h.sat_rich_SHEXI = LiBrSol.Calc_h_from_T_X_LiBrSol_Patek(T.sat_rich_SHEXI,x.LiBr_rich)
    if (h.sat_rich_SHEXI>h.sol_des_inI): # Saturation is not reached in SHEXI - no evaporation
        cp.sol_pump_out_molI = LiBrSol.Calc_cp_from_T_X_LiBrSol_Patek(T.sol_pump_outI,x.LiBr_rich)
        cp.sol_pump_outI = cp.sol_pump_out_molI/(x.LiBr_rich*M_LiBr + x.H2O_rich*M_H2O)
        T.sol_des_inI = T.sol_pump_outI + (h.sol_des_inI-h.sol_pump_outI)/cp.sol_pump_outI
    else:
        T.sol_des_inI = LiBrSol.Calc_state_SHEX_exit(T.sat_rich_SHEXI, p.cond_int, h.sol_des_inI, m.sol_richI, w.H2O_rich, 0.01)

    ## Poor solution after valve
    h.sol_abs_in = h.sol_valve_in
    T.sat_poor_valve = LiBrSol.Calc_T_from_p_X_satLiBrSol_Patek(p.evap,x.LiBr_poor)
    h.sat_poor_valve = LiBrSol.Calc_h_from_T_X_LiBrSol_Patek(T.sat_poor_valve,x.LiBr_poor)
    if (h.sat_poor_valve>h.sol_abs_in): # Saturation is not reached after valve - no evaporation
        T.sol_abs_in = T.sol_valve_in
    else:
        T.sol_abs_in = LiBrSol.Calc_state_valve_exit(T.sat_poor_valve, p.evap, h.sol_abs_in, m.sol_poor, w.H2O_poor, 0.01)
   
    ## Poor solution I after valve
    h.sol_abs_inI = h.sol_valve_inI
    T.sat_poor_valveI = LiBrSol.Calc_T_from_p_X_satLiBrSol_Patek(p.evap,x.LiBr_poorI)
    h.sat_poor_valveI = LiBrSol.Calc_h_from_T_X_LiBrSol_Patek(T.sat_poor_valveI,x.LiBr_poorI)
    if (h.sat_poor_valveI>h.sol_abs_inI): # Saturation is not reached after valve - no evaporation
        T.sol_abs_inI = T.sol_valve_inI
    else:
        T.sol_abs_inI = LiBrSol.Calc_state_valve_exit(T.sat_poor_valveI, p.evap, h.sol_abs_inI, m.sol_poorI, w.H2O_poorI, 0.01)

    #-------------------------------------------------------------------------#
    ## Post processing
    # Calculate fluxes over system boundary
    Q.cond = h.ref_cond_out*(m.ref+m.refI) - h.ref_cond_in*m.ref - h.ref_valve_outI*m.refI 
    Q.evap = (m.ref+m.refI)*(h.ref_evap_out-h.ref_evap_in)
    Q.abs = (m.sol_rich+m.sol_richI)*h.sol_abs_out - (m.ref+m.refI)*h.ref_abs_in - m.sol_poor*h.sol_abs_in - m.sol_poorI*h.sol_abs_inI
    Q.des =  m.refI*h.ref_des_outI + m.sol_poorI*h.sol_valve_inI - m.sol_richI*h.sol_pump_outI
    PP.W_pump = m.sol_rich*(h.sol_pump_out-h.sol_abs_out)
    PP.W_pumpI = m.sol_richI*(h.sol_pump_outI-h.sol_abs_out)
    # Heat Exchanger
    Q.RHEX = (m.ref+m.refI)*(h.ref_abs_in-h.ref_evap_out)
    h.RHEXideal = CoolProp.PropsSI('H','T',T.ref_cond_out,'P',p.evap,'Water')
    eta.RHEX = (h.ref_abs_in-h.ref_evap_out)/(h.RHEXideal-h.ref_evap_out)
    eta.SHEX = (h.sol_valve_in-h.sol_pump_out)/(h.sol_des_out-h.sol_pump_out)
    eta.SHEXI = (h.sol_valve_inI-h.sol_pump_outI)/(h.sol_des_outI-h.sol_pump_outI)
    # COP
    PP.COP = Q.evap/(Q.des + PP.W_pump + PP.W_pumpI)
    # Throttle loss
    PP.my_throttle = CoolProp.PropsSI('Q','H',h.ref_valve_in,'P',p.evap,'Water')
    # Circulation
    PP.f = (m.sol_rich+m.sol_richI)/(m.ref+m.refI)
    # Energy balance
    PP.energyBalance = Q.des + Q.des_mid + Q.evap + PP.W_pump + PP.W_pumpI + Q.cond + Q.condI + Q.abs
    # Mass balance (Absorber)
    PP.massBalance = (m.ref+m.refI) + (m.sol_poor+m.sol_poorI) - (m.sol_rich+m.sol_richI)
    # ----------------------------------------------------------------------- #
    ## Check
    # Refrigerant concentrations
    if (w.H2O_rich < w.H2O_poor):
        print("w_H2O_rich < w_H2O_poor")
        #error("w_H2O_rich < w_H2O_poor")

    if (w.H2O_rich < w.H2O_poorI):
        print("w_H2O_richI < w_H2O_poorI")
        #error("w_H2O_richI < w_H2O_poorI")

    if (w.H2O_rich < 0):
        print("w_H2O_rich < 0")
        #error("w_H2O_rich < 0")

    if (w.H2O_poor < 0):
        print("w_H2O_poor < 0")
        #error("w_H2O_poor < 0")

    if (w.H2O_poorI < 0):
        print("w_H2O_poorI < 0")
        #error("w_H2O_poorI < 0")

    if (w.H2O_rich - w.H2O_poor < 0.005):
        print("w_H2O_rich - w_H2O_poor < 0.005")
        #error("w_H2O_rich - w_H2O_poor < 0.005")

    if (w.H2O_rich - w.H2O_poorI < 0.005):
        print("w_H2O_rich - w_H2O_poorI < 0.005")
       #error("w_H2O_rich - w_H2O_poorI < 0.005")

    # Mass flow
    if (m.ref<0 or m.sol_poor<0 or m.sol_rich<0):
        print("mass flow is negativ")
        #error("mass flow is negativ")

    if (m.refI<0 or m.sol_poorI<0 or m.sol_richI<0):
        print("mass flow is negativ")
        #error("mass flow is negativ")

    # Crystallization and Violation where Patek is used
    # Weak Solution
    LiBrSol.checkForViolation_H2OLiBr(w.LiBr_rich,T.sol_abs_out,"Absorber exit")
    LiBrSol.checkForViolation_H2OLiBr(w.LiBr_rich,T.sol_pump_out,"Pump exit")
    LiBrSol.checkForViolation_H2OLiBr(w.LiBr_rich,T.sol_pump_outI,"PumpI exit")
    LiBrSol.checkForViolation_H2OLiBr(w.LiBr_rich,T.sol_des_in,"SHEX exit rich sol.")
    LiBrSol.checkForViolation_H2OLiBr(w.LiBr_rich,T.sol_des_inI,"SHEXI exit rich sol.")
    # Strong solution
    LiBrSol.checkForViolation_H2OLiBr(w.LiBr_poor,T.sol_des_out,"Desorber exit")
    LiBrSol.checkForViolation_H2OLiBr(w.LiBr_poorI,T.sol_des_outI,"DesorberI exit")
    LiBrSol.checkForViolation_H2OLiBr(w.LiBr_poor,T.sol_valve_in,"SHEX exit poor sol.")
    LiBrSol.checkForViolation_H2OLiBr(w.LiBr_poorI,T.sol_valve_inI,"SHEXI exit poor sol.")
    LiBrSol.checkForViolation_H2OLiBr(w.LiBr_poor,T.sol_abs_in,"Valve exit")
    LiBrSol.checkForViolation_H2OLiBr(w.LiBr_poorI,T.sol_abs_inI,"ValveI exit")
    # Energy and mass balance
    if (abs(PP.energyBalance) > 1):
        print('Energy is not conserved')
        #error("Energy is not conserved")

    if (abs(PP.massBalance) > 1):
        print('Mass is not conserved')
        #error("Mass is not conserved")

    #-------------------------------------------------------------------------#

    return T, p, h, m, w, eta, Q, PP, s