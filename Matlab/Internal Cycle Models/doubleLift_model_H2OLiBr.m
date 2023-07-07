function [T, p, h, m, w, eta, Q, PP, s] = doubleLift_model_H2OLiBr(T, p, h, m, eta, Q, HX, s)
%% Function Double Lift Model H2OLiBr
% ----------------------------------------------------------------------- %
%{
Author  : Ludwig Irrgang
Date    : 01.07.2023
Copyright information:
Ludwig Irrgang
Lehrstuhl für Energiesysteme
TUM School of Engineering and Design
Technische Universität München
Boltzmannstr. 15 
85748 Garching b. München
ludwig.irrgang@tum.de
%}
% ----------------------------------------------------------------------- %
if nargin<8||isempty(s),error('Input Argument: Setup s missing');end
if nargin<7||isempty(HX),error('Input Argument: Approach temperature missing');end
if nargin<6||isempty(Q),error('Input Argument: Heat missing');end
if nargin<5||isempty(eta),error('Input Argument: Efficiency missing');end
if nargin<4||isempty(m),error('Input Argument: Mass missing');end
if nargin<3||isempty(h),error('Input Argument: Enthalpie missing');end
if nargin<2||isempty(p),error('Input Argument: Pressure missing');end
if nargin<1||isempty(T),error('Input Argument: Temperature missing');end
% ----------------------------------------------------------------------- %
%% Input/ Output
% ----------------------------------------------------------------------- %
% Input:
%{
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
%}
% Output:
%{
1. Temperature Struct                       [K]
2. Pressure Struct                          [Pa]
3. Specific enthalpy Struct                 [J/kg]
4. Mass flow rate Struct                    [kg/s]
5. Mass fraction Struct                     [kg/s]
6. Efficiency Struct                        [-]
7. Heat flow Struct                         [J]
8. Post Process Struct
9. Setup struct
%}
% ----------------------------------------------------------------------- %
%% Absorption System
% ----------------------------------------------------------------------- %
% Components:
%{
- 2 x Desorber
- 2 x Absorber
- 2 x Pump
- 3 x Throttle valves
- Condenser
- Evaporator
- Solution heat exchanger
- Refrigerant heat exchanger
- Working fluid: Water lithium-bromide solution
- Refrigerant: Water
%}
% ----------------------------------------------------------------------- %
%% Annotations
% ----------------------------------------------------------------------- %
%{ 
- x stands for molar fraction
- w stands for mass fraction
- Functions for LiBr are using molar state properties and molar fraction
%} 
% ----------------------------------------------------------------------- %
%% Assumptions
% ----------------------------------------------------------------------- %
%{ 
- Temperature of refrigerant leaving desorber is 5K below desorber temp.
- All components of the system operate in steady state
- Solution leaving generator and absorber is saturated
- Refrigerant leaving evaporator is saturated
- Pressure drops in the system components are negelcted
- Heat capacity of solution assumed to be constant in undersaturated solutions
%} 
% ----------------------------------------------------------------------- %
%% Define of constants
% ------------------------Necessary Constants---------------------------- %
M_LiBr = 0.08685;                 %[kg/mol]
M_H2O = 0.018015268;              %[kg/mol]
% ----------------------------------------------------------------------- %
%% Calculation
% ---------------------------CALCULATION--------------------------------- %
%% Calculate pressures
p.evap = CoolProp.PropsSI('P','T',T.evap,'Q',1,'Water');
p.cond = CoolProp.PropsSI('P','T',T.cond,'Q',0,'Water');
p.mid = sqrt(p.evap*p.cond);
% ----------------------------------------------------------------------- %
%% Rich solution (High ref. concentration)
% Low pressure cycle
% Absorber
x.LiBr_rich = Calc_X_from_T_p_satLiBrSol_Patek(T.sol_abs_out,p.evap);
x.H2O_rich = 1 - x.LiBr_rich;
w.LiBr_rich = x.LiBr_rich*M_LiBr/ (x.LiBr_rich*M_LiBr + x.H2O_rich*M_H2O);
w.H2O_rich = 1 - w.LiBr_rich;
h.sol_abs_out_mol = Calc_h_from_T_X_LiBrSol_Patek(T.sol_abs_out,x.LiBr_rich);
h.sol_abs_out = h.sol_abs_out_mol/(x.LiBr_rich*M_LiBr + x.H2O_rich*M_H2O);
rho.sol_abs_out_mol = Calc_rho_from_T_X_LiBrSol_Patek(T.sol_abs_out,x.LiBr_rich);
rho.sol_abs_out = rho.sol_abs_out_mol*(x.LiBr_rich*M_LiBr + x.H2O_rich*M_H2O);
v.sol_abs_out = 1/rho.sol_abs_out;
% Pump
w.pump = v.sol_abs_out*(p.mid - p.evap);
h.sol_pump_out = h.sol_abs_out + w.pump/eta.pump;
% ----------------------------------------------------------------------- %
% Prohibit crystallization in SHEXI - necessary internal condensation pressure increase 
T.sol_pump_outI = T.sol_abs_outI; % Is assumed to be isothermal
% SHEX
if(T.sol_pump_outI + HX.T_PP_SHEXI < T.sol_des_outI)
    T.sol_valve_inI = T.sol_pump_outI + HX.T_PP_SHEXI;
else
    T.sol_valve_inI = T.sol_des_outI;
end
w.crI = crystallization_H2OLiBr("T",T.sol_valve_inI);
x.crI = w.crI/M_LiBr/ (w.crI/M_LiBr+(1-w.crI)/M_H2O);
p.crI = Calc_p_from_T_X_LiBrSol_Patek(T.sol_des_outI,x.crI);
if (p.cond<p.crI)
    p.cond = p.crI + 10; % Prevent calculation error
    T.cond = CoolProp.PropsSI('T','Q',0,'P',p.cond,'Water');
end
% ----------------------------------------------------------------------- %
% Prohibit crystallization in SHEX - necessary desorber temperature
% decrease, since pressure is already fixed for maximal efficiency
T.sol_pump_out = T.sol_abs_out; % Is assumed to be isothermal
% SHEX
if(T.sol_pump_out + HX.T_PP_SHEX < T.sol_des_out)
    T.sol_valve_in = T.sol_pump_out + HX.T_PP_SHEX;
else
    T.sol_valve_in = T.sol_des_out;
end
w.cr = crystallization_H2OLiBr("T",T.sol_valve_in);
x.cr = w.cr/M_LiBr/ (w.cr/M_LiBr+(1-w.cr)/M_H2O);
p.cr = Calc_p_from_T_X_LiBrSol_Patek(T.sol_des_out,x.cr);
while(p.mid<p.cr)
    if(T.sol_pump_out + HX.T_PP_SHEX < T.sol_des_out)
        T.sol_valve_in = T.sol_pump_out + HX.T_PP_SHEX;
    else
        T.sol_valve_in = T.sol_des_out;
    end
    w.cr = crystallization_H2OLiBr("T",T.sol_valve_in);
    x.cr = w.cr/M_LiBr/ (w.cr/M_LiBr+(1-w.cr)/M_H2O);
    p.cr = Calc_p_from_T_X_LiBrSol_Patek(T.sol_des_out,x.cr);
    if (p.mid<p.cr)
        T.sol_des_out = T.sol_des_out - 0.01; % Iterate until temperature low enough
    end
end
% ----------------------------------------------------------------------- %
% High pressure cycle
% Absorber
x.LiBr_richI = Calc_X_from_T_p_satLiBrSol_Patek(T.sol_abs_outI,p.mid);
x.H2O_richI = 1 - x.LiBr_richI;
w.LiBr_richI = x.LiBr_richI*M_LiBr/ (x.LiBr_richI*M_LiBr + x.H2O_richI*M_H2O);
w.H2O_richI = 1 - w.LiBr_richI;
h.sol_abs_out_molI = Calc_h_from_T_X_LiBrSol_Patek(T.sol_abs_outI,x.LiBr_richI);
h.sol_abs_outI = h.sol_abs_out_molI/(x.LiBr_richI*M_LiBr + x.H2O_richI*M_H2O);
rho.sol_abs_out_molI = Calc_rho_from_T_X_LiBrSol_Patek(T.sol_abs_outI,x.LiBr_richI);
rho.sol_abs_outI = rho.sol_abs_out_molI*(x.LiBr_richI*M_LiBr + x.H2O_richI*M_H2O);
v.sol_abs_outI = 1/rho.sol_abs_outI;
% Pump
w.pumpI = v.sol_abs_outI*(p.cond - p.mid);
h.sol_pump_outI = h.sol_abs_outI + w.pumpI/eta.pump;
%-------------------------------------------------------------------------%
%% Refrigerant line
% Desorber
T.ref_des_out = T.sol_des_out - HX.dT_ref_des;
h.ref_des_out = CoolProp.PropsSI('H','T',T.ref_des_out,'P',p.mid,'Water');
T.ref_des_outI = T.sol_des_outI - HX.dT_ref_desI;
h.ref_des_outI = CoolProp.PropsSI('H','T',T.ref_des_outI,'P',p.cond,'Water');
% Condenser
T.ref_cond_in = T.ref_des_outI;
h.ref_cond_in = h.ref_des_outI;
if (T.ext_cond_in+HX.T_PP_cond<T.cond)
    T.ref_cond_out = T.ext_cond_in + HX.T_PP_cond; % Subcooling as low as possible
    h.ref_cond_out = CoolProp.PropsSI('H','P',p.cond,'T',T.ref_cond_out,'Water');
else
    T.ref_cond_out = T.cond;
    h.ref_cond_out = CoolProp.PropsSI('H','P',p.cond,'Q',0,'Water');
end
% Evaporator
T.ref_evap_out = T.evap;
h.ref_evap_out = CoolProp.PropsSI('H','P',p.evap,'Q',1,'Water');
% Subcooler (heat capacity of steam lower than liquid)
if(T.ref_cond_out - HX.T_PP_RHEX > T.ref_evap_out)
    T.ref_abs_in = T.ref_cond_out - HX.T_PP_RHEX;
    h.ref_abs_in = CoolProp.PropsSI('H','T',T.ref_abs_in,'P',p.evap,'Water'); 
else
    T.ref_abs_in = T.ref_evap_out;
    h.ref_abs_in = h.ref_evap_out; 
end
h.ref_valve_in = h.ref_cond_out - (h.ref_abs_in-h.ref_evap_out);
T.ref_valve_in = CoolProp.PropsSI('T','H',h.ref_valve_in,'P',p.cond,'Water');
% Expansion valve
h.ref_evap_in = h.ref_valve_in; % Isenthalpic thottle
T.ref_evap_in = CoolProp.PropsSI('T','H',h.ref_evap_in,'P',p.evap,'Water');
% ----------------------------------------------------------------------- %
%% Poor solution (Low ref. concentration)
% Low pressure
% Desorber
x.LiBr_poor = Calc_X_from_T_p_satLiBrSol_Patek(T.sol_des_out,p.mid);
x.H2O_poor = 1 - x.LiBr_poor;
w.LiBr_poor = x.LiBr_poor*M_LiBr/ (x.LiBr_poor*M_LiBr + x.H2O_poor*M_H2O);
w.H2O_poor = 1 - w.LiBr_poor;
h.sol_des_out_mol = Calc_h_from_T_X_LiBrSol_Patek(T.sol_des_out,x.LiBr_poor);
h.sol_des_out = h.sol_des_out_mol/(x.LiBr_poor*M_LiBr + x.H2O_poor*M_H2O);
% Assuming constant heat capacity
cp.sol_des_out_mol = Calc_cp_from_T_X_LiBrSol_Patek(T.sol_des_out,x.LiBr_poor);
cp.sol_des_out = cp.sol_des_out_mol/(x.LiBr_poor*M_LiBr + x.H2O_poor*M_H2O);
h.sol_valve_in = h.sol_des_out+cp.sol_des_out*(T.sol_valve_in-T.sol_des_out);
% High pressure
% Desorber
x.LiBr_poorI = Calc_X_from_T_p_satLiBrSol_Patek(T.sol_des_outI,p.cond);
x.H2O_poorI = 1 - x.LiBr_poorI;
w.LiBr_poorI = x.LiBr_poorI*M_LiBr/ (x.LiBr_poorI*M_LiBr + x.H2O_poorI*M_H2O);
w.H2O_poorI = 1 - w.LiBr_poorI;
h.sol_des_out_molI = Calc_h_from_T_X_LiBrSol_Patek(T.sol_des_outI,x.LiBr_poorI);
h.sol_des_outI = h.sol_des_out_molI/(x.LiBr_poorI*M_LiBr + x.H2O_poorI*M_H2O);
% Assuming constant heat capacity
cp.sol_des_out_molI = Calc_cp_from_T_X_LiBrSol_Patek(T.sol_des_outI,x.LiBr_poorI);
cp.sol_des_outI = cp.sol_des_out_molI/(x.LiBr_poorI*M_LiBr + x.H2O_poorI*M_H2O);
h.sol_valve_inI = h.sol_des_outI+cp.sol_des_outI*(T.sol_valve_inI-T.sol_des_outI);
% ----------------------------------------------------------------------- %
%% Energy balances and mass conservation
% Solve system of linear equations (energy conservation, mass conservation)
% Solution vector: Q.des, m.ref_des_out, m.sol_rich, m.sol_poor, Q.desI,
% m.refI, m.sol_richI, m.sol_poorI
switch s.requirement
    case 'Q_des'
        A = [   1,  0,              0,              0,                  1,          0,                  0,                  0;
                0,  -1,             1,              -1,                 0,          0,                  0,                  0;
                1,  -h.ref_des_out, h.sol_pump_out, -h.sol_valve_in,    0,          0,                  0,                  0;
                0,  -1,             w.H2O_rich,     -w.H2O_poor,        0,          0,                  0,                  0;
                0,  0,              0,              0,                  0,          -1,                 1,                  -1;
                0,  0,              0,              0,                  1,          -h.ref_des_outI,    h.sol_pump_outI,    -h.sol_valve_inI;
                0,  0,              0,              0,                  0,          -1,                 w.H2O_richI,        -w.H2O_poorI;
                0,  1,              0,              0,                  0,          -1,                 0,                  0   ];
        b = [   Q.dec;  0;  0;  0;  0;  0;  0;  0  ];
    case 'Q_evap'
    A = [       0,  0,              0,              0,                  0,          (h.ref_evap_out-h.ref_evap_in), 0,                  0;
                0,  -1,             1,              -1,                 0,          0,                              0,                  0;
                1,  -h.ref_des_out, h.sol_pump_out, -h.sol_valve_in,    0,          0,                              0,                  0;
                0,  -1,             w.H2O_rich,     -w.H2O_poor,        0,          0,                              0,                  0;
                0,  0,              0,              0,                  0,          -1,                             1,                  -1;
                0,  0,              0,              0,                  1,          -h.ref_des_outI,                h.sol_pump_outI,    -h.sol_valve_inI;
                0,  0,              0,              0,                  0,          -1,                             w.H2O_richI,        -w.H2O_poorI;
                0,  1,              0,              0,                  0,          -1,                             0,                  0   ];
        b = [   Q.dec;  0;  0;  0;  0;  0;  0;  0  ];
    otherwise
        error('AKM decider is not defined properly. Use Q_des or Q_evap')
end
y = A\b;
% Low pressure side
Q.des = y(1);
m.ref = y(2);
m.sol_rich = y(3);
m.sol_poor = y(4);
% High pressure side
Q.desI = y(5);
m.refI = y(6);
m.sol_richI = y(7);
m.sol_poorI = y(8);
% ----------------------------------------------------------------------- %
%% Rich solution after SHEX
Q.SHEX = (h.sol_des_out - h.sol_valve_in)*m.sol_poor;
h.sol_des_in = h.sol_pump_out + Q.SHEX/m.sol_rich;
cp.sol_pump_out_mol = Calc_cp_from_T_X_LiBrSol_Patek(T.sol_pump_out,x.LiBr_rich);
cp.sol_pump_out = cp.sol_pump_out_mol/(x.LiBr_rich*M_LiBr + x.H2O_rich*M_H2O);
T.sol_des_in = T.sol_pump_out + (h.sol_des_in-h.sol_pump_out)/cp.sol_pump_out;
%% Rich solution after SHEXI
Q.SHEXI = (h.sol_des_outI - h.sol_valve_inI)*m.sol_poorI;
h.sol_des_inI = h.sol_pump_outI + Q.SHEXI/m.sol_richI;
cp.sol_pump_out_molI = Calc_cp_from_T_X_LiBrSol_Patek(T.sol_pump_outI,x.LiBr_richI);
cp.sol_pump_outI = cp.sol_pump_out_molI/(x.LiBr_richI*M_LiBr + x.H2O_richI*M_H2O);
T.sol_des_inI = T.sol_pump_outI + (h.sol_des_inI-h.sol_pump_outI)/cp.sol_pump_outI;
%% Poor solution after valve
h.sol_abs_in = h.sol_valve_in;
T.sol_abs_in = T.sol_valve_in;
%% Poor solution I after valve
h.sol_abs_inI = h.sol_valve_inI;
T.sol_abs_inI = T.sol_valve_inI;
% ----------------------------------------------------------------------- %
%% Post processing
% Fluxes over system boundary
Q.des_ges = Q.des + Q.desI;
Q.cond = (h.ref_cond_out-h.ref_cond_in)*m.refI;
Q.evap = (h.ref_evap_out-h.ref_evap_in)*m.refI;
Q.abs =  h.sol_abs_out*m.sol_rich - h.sol_abs_in*m.sol_poor - h.ref_abs_in*m.refI;
Q.absI = h.sol_abs_outI*m.sol_richI - h.sol_abs_inI*m.sol_poorI - h.ref_des_out*m.ref;
PP.W_pump = m.sol_rich*(h.sol_pump_out - h.sol_abs_out);
PP.W_pumpI = m.sol_richI*(h.sol_pump_outI - h.sol_abs_outI);
% Heat Exchanger
Q.RHEX = m.refI*(h.ref_abs_in-h.ref_evap_out);
h.RHEXideal = CoolProp.PropsSI('H','T',T.ref_cond_out,'P',p.evap,'Water');
eta.RHEX = (h.ref_abs_in-h.ref_evap_out)/(h.RHEXideal-h.ref_evap_out);
eta.SHEX = (h.sol_valve_in-h.sol_pump_out)/(h.sol_des_out-h.sol_pump_out);
eta.SHEXI = (h.sol_valve_inI-h.sol_pump_outI)/(h.sol_des_outI-h.sol_pump_outI);
% COP
PP.COP = Q.evap/(Q.des_ges + PP.W_pump + PP.W_pumpI);
% Throttle loss
PP.my_throttle = CoolProp.PropsSI('Q','H',h.ref_valve_in,'P',p.evap,'Water');
% Circulation
PP.f = (m.sol_rich+m.sol_richI)/(m.ref+m.refI);
% Energy balance
PP.energyBalance = Q.des_ges + Q.evap + PP.W_pump + PP.W_pumpI + Q.cond + Q.abs + Q.absI;
% Mass balance
PP.massBalance = m.sol_rich - m.sol_poor - m.refI;
% ----------------------------------------------------------------------- %
%% Check
% Refrigerant concentrations
if (w.H2O_rich < w.H2O_poor)
    error("w_H2O_rich < w_H2O_poor")
end
if (w.H2O_richI < w.H2O_poorI)
    error("w_H2O_richI < w_H2O_poorI")
end
if (w.H2O_rich < 0)
    error("w_H2O_rich < 0")
end
if (w.H2O_richI < 0)
    error("w_H2O_richI < 0")
end
if (w.H2O_poor < 0)
    error("w_H2O_poor < 0")
end
if (w.H2O_poorI < 0)
    error("w_H2O_poorI < 0")
end
% Mass flow
if (m.ref<0 || m.sol_poor<0 || m.sol_rich<0)
        error("mass flow is negativ")
end
if (m.refI<0 || m.sol_poorI<0 || m.sol_richI<0)
        error("mass flow is negativ")
end
% Crystallization and Violation where Patek is used
% Rich Solution
checkForViolation_H2OLiBr(w.LiBr_rich,T.sol_abs_out,"Absorber exit");
checkForViolation_H2OLiBr(w.LiBr_richI,T.sol_abs_outI,"Absorber I exit");
checkForViolation_H2OLiBr(w.LiBr_rich,T.sol_pump_out,"Pump exit");
checkForViolation_H2OLiBr(w.LiBr_richI,T.sol_pump_outI,"PumpI exit");
checkForViolation_H2OLiBr(w.LiBr_rich,T.sol_des_in,"SHEX exit rich sol.");
checkForViolation_H2OLiBr(w.LiBr_richI,T.sol_des_inI,"SHEXI exit rich sol.");
% Poor solution
checkForViolation_H2OLiBr(w.LiBr_poor,T.sol_des_out,"Desorber exit");
checkForViolation_H2OLiBr(w.LiBr_poorI,T.sol_des_outI,"DesorberI exit");
checkForViolation_H2OLiBr(w.LiBr_poor,T.sol_valve_in,"SHEX exit poor sol.");
checkForViolation_H2OLiBr(w.LiBr_poorI,T.sol_valve_inI,"SHEXI exit poor sol.");
checkForViolation_H2OLiBr(w.LiBr_poor,T.sol_abs_in,"Valve exit");
checkForViolation_H2OLiBr(w.LiBr_poorI,T.sol_abs_inI,"ValveI exit");
% Energy and mass balance
if (abs(PP.energyBalance) > 10)
    error("Energy is not conserved")
end
if (abs(PP.massBalance) > 1)
    error("Mass is not conserved")
end
% ----------------------------------------------------------------------- %
end
