function [T, p, h, m, w, eta, Q, PP, s] = base_model_H2OLiBr(T, p, h, m, eta, Q, HX, s)
%% Function Base Model H20LiBr
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
- HX.T_PP_RHEX
- HX.T_PP_cond
- HX.dT_ref_des
%}
% Output:
%{
1. Temperature struct                       [K]
2. Pressure struct                          [Pa]
3. Specific enthalpy struct                 [J/kg]
4. Mass flow rate struct                    [kg/s]
5. Mass fraction struct                     [kg/s]
6. Efficiency struct                        [-]
7. Heat flow struct                         [J]
8. Post process struct
9. Setup/ Entropy struct
%}
% ----------------------------------------------------------------------- %
%% Absorption System
% ----------------------------------------------------------------------- %
% Components:
%{
- Desorber
- Condenser
- Evaporator
- Absorber
- Pump
- 2 x Throttle Valves
- Solution heat exchanger
- Refrigerant heat exchanger
- Working fluid: water lithium-bromide solution
- Refrigerant: water
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
%% Definition of constants
% ------------------------Necessary Constants---------------------------- %
M_LiBr = 0.08685;               %[kg/mol]
M_H2O = 0.018015268;            %[kg/mol]
%% Calculation
% ---------------------------CALCULATION--------------------------------- %
%% Calculate pressures
p.evap = CoolProp.PropsSI('P','T',T.evap,'Q',1,'Water');
p.cond = CoolProp.PropsSI('P','T',T.cond,'Q',0,'Water');
%-------------------------------------------------------------------------%
%% Rich solution (High ref. concentration)
% Absorber
x.LiBr_rich = Calc_X_from_T_p_satLiBrSol_Patek(T.sol_abs_out,p.evap);
x.H2O_rich = 1 - x.LiBr_rich;
w.LiBr_rich = x.LiBr_rich*M_LiBr/ (x.LiBr_rich*M_LiBr + x.H2O_rich*M_H2O);
w.H2O_rich = 1 - w.LiBr_rich;
h.sol_abs_out_mol = Calc_h_from_T_X_LiBrSol_Patek(T.sol_abs_out,x.LiBr_rich);
h.sol_abs_out = h.sol_abs_out_mol/(x.LiBr_rich*M_LiBr + x.H2O_rich*M_H2O);
rho.sol_abs_out_mol = Calc_rho_from_T_X_LiBrSol_Patek(T.sol_abs_out,x.LiBr_rich);
rho.sol_abs_out = rho.sol_abs_out_mol*(x.LiBr_rich*M_LiBr + x.H2O_rich*M_H2O);
% ----------------------------------------------------------------------- %
% Prohibit crystallization in SHEX - necessary condensation pressure increase 
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
if (p.cond<p.cr)
    p.cond = p.cr + 10; % Prevent calculation error and add safety 
    T.cond = CoolProp.PropsSI('T','Q',0,'P',p.cond,'Water');
end
% ----------------------------------------------------------------------- %
% Pump
h.sol_pump_out = h.sol_abs_out + (1/rho.sol_abs_out)*(p.cond - p.evap)/eta.pump;
%-------------------------------------------------------------------------%
%% Refrigerant line
% Desorber
T.ref_des_out = T.sol_des_out - HX.dT_ref_des;
h.ref_des_out = CoolProp.PropsSI('H','T',T.ref_des_out,'P',p.cond,'Water');
% Condenser
T.ref_cond_in = T.ref_des_out;
h.ref_cond_in = h.ref_des_out;
if (T.ext_cond_in+HX.T_PP_cond<T.cond)
    T.ref_cond_out = T.ext_cond_in + HX.T_PP_cond; % Sub-cooling as low as possible
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
% Throttle
h.ref_evap_in = h.ref_valve_in; % Isenthalpic thottle
T.ref_evap_in = CoolProp.PropsSI('T','H',h.ref_evap_in,'P',p.evap,'Water');
%-------------------------------------------------------------------------%
%% Poor solution (Low ref. concentration)
% Desorber
x.LiBr_poor = Calc_X_from_T_p_satLiBrSol_Patek(T.sol_des_out,p.cond);
x.H2O_poor = 1 - x.LiBr_poor;
w.LiBr_poor = x.LiBr_poor*M_LiBr/ (x.LiBr_poor*M_LiBr + x.H2O_poor*M_H2O);
w.H2O_poor = 1 - w.LiBr_poor;
h.sol_des_out_mol = Calc_h_from_T_X_LiBrSol_Patek(T.sol_des_out,x.LiBr_poor);
h.sol_des_out = h.sol_des_out_mol/(x.LiBr_poor*M_LiBr + x.H2O_poor*M_H2O);
% Assuming constant heat capacity
cp.sol_des_out_mol = Calc_cp_from_T_X_LiBrSol_Patek(T.sol_des_out,x.LiBr_poor);
cp.sol_des_out = cp.sol_des_out_mol/(x.LiBr_poor*M_LiBr + x.H2O_poor*M_H2O);
h.sol_valve_in = h.sol_des_out+cp.sol_des_out*(T.sol_valve_in-T.sol_des_out);
%-------------------------------------------------------------------------%
%% Energy balances and mass conservation
% Solve system of linear equations (energy conservation, mass conservation)
% Solution vector: m_sol_rich, m_ref, m_sol_poor
switch s.requirement
    case 'Q_des'
        A = [   h.sol_pump_out,     -h.ref_des_out,     -h.sol_valve_in;
                w.H2O_rich,         -1,                 -w.H2O_poor;
                1,                  -1,                 -1   ];
        b = [   -Q.dec;             0;                  0    ];
    case 'Q_evap'
        A = [   0,          (h.ref_evap_out-h.ref_evap_in),     0;
                w.H2O_rich, -1,                                 -w.H2O_poor;
                1,          -1,                                 -1   ];
        b = [   Q.dec;      0;                                  0    ];
    otherwise
        error('AKM requirement is not defined properly. Use Q_des or Q_evap')
end
y = A\b;
m.sol_rich = y(1);
m.ref = y(2);
m.sol_poor = y(3);
%-------------------------------------------------------------------------%
%% Rich solution after SHEX
Q.SHEX = (h.sol_des_out - h.sol_valve_in)*m.sol_poor;
h.sol_des_in = h.sol_pump_out + Q.SHEX/m.sol_rich;
cp.sol_pump_out_mol = Calc_cp_from_T_X_LiBrSol_Patek(T.sol_pump_out,x.LiBr_rich);
cp.sol_pump_out = cp.sol_pump_out_mol/(x.LiBr_rich*M_LiBr + x.H2O_rich*M_H2O);
% Assuming constant heat capacity
T.sol_des_in = T.sol_pump_out + (h.sol_des_in-h.sol_pump_out)/cp.sol_pump_out;
%% Poor solution after valve
h.sol_abs_in = h.sol_valve_in;
T.sol_abs_in = T.sol_valve_in;
%-------------------------------------------------------------------------%
%% Post processing
% Fluxes over system boundary
Q.cond = m.ref*(h.ref_cond_out-h.ref_cond_in);
Q.evap = m.ref*(h.ref_evap_out-h.ref_evap_in);
Q.abs =  m.sol_rich*h.sol_abs_out - m.ref*h.ref_abs_in - m.sol_poor*h.sol_abs_in;
PP.W_pump = m.sol_rich*(h.sol_pump_out-h.sol_abs_out);
Q.des = m.ref*h.ref_des_out + m.sol_poor*h.sol_valve_in - m.sol_rich*h.sol_des_in;
% Heat Exchanger
Q.SHEX = m.sol_poor*(h.sol_des_out-h.sol_valve_in);
Q.RHEX = m.ref*(h.ref_abs_in-h.ref_evap_out);
h.RHEXideal = CoolProp.PropsSI('H','T',T.ref_cond_out,'P',p.evap,'Water');
eta.RHEX = (h.ref_abs_in-h.ref_evap_out)/(h.RHEXideal-h.ref_evap_out);
eta.SHEX = (h.sol_des_in-h.sol_pump_out)/(h.sol_des_out-h.sol_pump_out);
% COP
PP.COP = Q.evap/(Q.des + PP.W_pump);
% Throttle loss
PP.my_throttle = CoolProp.PropsSI('Q','H',h.ref_valve_in,'P',p.evap,'Water');
% Circulation
PP.f = m.sol_rich/m.ref;
% Energy balance
PP.energyBalance = Q.des + Q.evap + PP.W_pump + Q.cond + Q.abs;
% Mass balance
PP.massBalance = m.ref + m.sol_poor - m.sol_rich;
%-------------------------------------------------------------------------%
%% Check
% Refrigerant concentrations
if (w.H2O_rich < w.H2O_poor)
    error("w_H2O_rich < w_H2O_poor")
end
if (w.H2O_rich < 0)
    error("w_H2O_rich < 0")
end
if (w.H2O_poor < 0)
    error("w_H2O_poor < 0")
end
if (w.H2O_rich - w.H2O_poor < 0.005)
    error("w_H2O_rich - w_H2O_poor < 0.005")
end
% Mass flow
if (m.ref<0 || m.sol_poor<0 || m.sol_rich<0)
        error("mass flow is negativ")
end
% Crystallization and Violation where Patek is used
% Rich Solution
checkForViolation_H2OLiBr(w.LiBr_rich,T.sol_abs_out,"Absorber exit");
checkForViolation_H2OLiBr(w.LiBr_rich,T.sol_pump_out,"Pump exit");
checkForViolation_H2OLiBr(w.LiBr_rich,T.sol_des_in,"SHEX exit rich sol.");
% Poor solution
checkForViolation_H2OLiBr(w.LiBr_poor,T.sol_des_out,"Desorber exit");
checkForViolation_H2OLiBr(w.LiBr_poor,T.sol_valve_in,"SHEX exit poor sol.");
checkForViolation_H2OLiBr(w.LiBr_poor,T.sol_abs_in,"Valve exit");
% Energy and mass balance
if (abs(PP.energyBalance) > 0.1)
    error("Energy is not conserved")
end
if (abs(PP.massBalance) > 1)
    error("Mass is not conserved")
end
%-------------------------------------------------------------------------%
end
