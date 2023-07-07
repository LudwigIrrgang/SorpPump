function [T, p, h, m, w, eta, Q, PP, s] = base_model_NH3H2O(T, p, h, m, eta, Q, HX, s)
%% Function Base Model NH3H2O
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
- Solution Heat Exchanger
- Refrigerant Heat Exchanger
- Working fluid: Ammonia Water Solution
- Refrigerant: Ammonia
%}
% ----------------------------------------------------------------------- %
%% Assumptions
% ----------------------------------------------------------------------- %
%{
- Temperature of refrigerant leaving desorber is 5K below desorber temp.
- All components of the system operate in steady state
- Solution leaving generator and absorber is saturated
- Refrigerant leaving condenser and evaporator is saturated
- Pressure drops in the system components are negelcted
%} 
% ----------------------------------------------------------------------- %
%% Calculation
% ---------------------------CALCULATION--------------------------------- %
%% Calculate pressures
% Low pressure
p.evap = CoolProp.PropsSI('P','T',T.evap,'Q',1,'AMMONIA');
p.cond = CoolProp.PropsSI('P','T',T.cond,'Q',0,'AMMONIA');
%% Refrigerant line
% Desorber
T.ref_des_out = T.sol_des_out - HX.dT_ref_des;
h.ref_des_out = CoolProp.PropsSI('H','T',T.ref_des_out,'P',p.cond,'AMMONIA');
% Condenser
T.ref_cond_in = T.ref_des_out;
h.ref_cond_in = h.ref_des_out;
if (T.ext_cond_in+HX.T_PP_cond<T.cond)
    T.ref_cond_out = T.ext_cond_in + HX.T_PP_cond; % Subcooling as low as possible
    h.ref_cond_out = CoolProp.PropsSI('H','T',T.ref_cond_out,'P',p.cond,'AMMONIA');
else
    T.ref_cond_out = T.cond;
    h.ref_cond_out = CoolProp.PropsSI('H','T',T.ref_cond_out,'Q',0,'AMMONIA');
end
% Evaporator
T.ref_evap_out = T.evap;
h.ref_evap_out = CoolProp.PropsSI('H','P',p.evap,'Q',1,'AMMONIA');
% Subcooler (heat capacity of steam lower than liquid)
if(T.ref_cond_out - HX.T_PP_RHEX > T.ref_evap_out)
    T.ref_abs_in = T.ref_cond_out - HX.T_PP_RHEX;
    h.ref_abs_in = CoolProp.PropsSI('H','T',T.ref_abs_in,'P',p.evap,'AMMONIA');
else
    T.ref_abs_in = T.ref_evap_out;
    h.ref_abs_in = h.ref_evap_out;
end
h.ref_valve_in = h.ref_cond_out - (h.ref_abs_in-h.ref_evap_out);
T.ref_valve_in = CoolProp.PropsSI('T','H',h.ref_valve_in,'P',p.cond,'AMMONIA');
% Throttle
h.ref_evap_in = h.ref_valve_in; % Isenthalpic throttle
T.ref_evap_in = CoolProp.PropsSI('T','H',h.ref_evap_in,'P',p.evap,'AMMONIA');
%-------------------------------------------------------------------------%
%% Rich solution (High ref. concentration)
% Absorber
w.NH3_rich = NH3inSolution_Calc_X_PT(p.evap/1000,T.sol_abs_out);
h.sol_abs_out = refpropm('H','T',T.sol_abs_out,'Q',0,'AMMONIA','WATER',[w.NH3_rich (1-w.NH3_rich)]);
% Pump
s.sol_abs_out = refpropm('S','T',T.sol_abs_out,'Q',0,'AMMONIA','WATER',[w.NH3_rich (1-w.NH3_rich)]);
s.sol_pump_out = s.sol_abs_out;
roh_sol_abs_out = refpropm('D','T',T.sol_abs_out,'Q',0,'AMMONIA','WATER',[w.NH3_rich (1-w.NH3_rich)]);
v_sol_abs_out = 1/roh_sol_abs_out;
h.sol_pump_out = h.sol_abs_out+v_sol_abs_out*(p.cond-p.evap)/eta.pump;
T.sol_pump_out = T.sol_abs_out;
cp.sol_pump_out = refpropm('C','T',T.sol_pump_out,'P',p.cond/1000,'AMMONIA','WATER',[w.NH3_rich (1-w.NH3_rich)]);
%-------------------------------------------------------------------------%
%% Poor solution (Low ref. concentration)
% Desorber
w.NH3_poor = NH3inSolution_Calc_X_PT(p.cond/1000,T.sol_des_out);
if(w.NH3_poor>0)
    h.sol_des_out = refpropm('H','T',T.sol_des_out,'Q',0,'AMMONIA','WATER',[w.NH3_poor (1-w.NH3_poor)]);
else
    h.sol_des_out = refpropm('H','T',T.sol_des_out,'P',p.cond,'WATER');
end
% SHEX
if(T.sol_pump_out + HX.T_PP_SHEX < T.sol_des_out)
    T.sol_valve_in = T.sol_pump_out + HX.T_PP_SHEX;
else
    T.sol_valve_in = T.sol_des_out;
end
h.sol_valve_in = refpropm('H','T',T.sol_valve_in,'P',p.cond/1000,'AMMONIA','WATER',[w.NH3_poor (1-w.NH3_poor)]);
%-------------------------------------------------------------------------%
%% Energy balances and mass conservation
% Solve system of linear equations (energy conservation, mass conservation)
% Solution vector: m_rich, m_ref, m_poor
switch s.requirement
    case 'Q_des'
        A = [   h.sol_pump_out,     -h.ref_des_out,     -h.sol_valve_in;
                w.NH3_rich,         -1,                 -w.NH3_poor;
                1,                  -1,                 -1   ];
        b = [   -Q.dec;     0;      0   ];
    case 'Q_evap'
        A = [   0,          (h.ref_evap_out-h.ref_evap_in), 0;
                w.NH3_rich, -1,                             -w.NH3_poor;
                1,          -1,                             -1   ];
        b = [   Q.dec;     0;      0   ];
    otherwise
        error('AKM requirement is not defined properly. Use Q_des or Q_evap')
end
x = A\b;
m.sol_rich = x(1);
m.ref = x(2);
m.sol_poor = x(3);
%-------------------------------------------------------------------------%
%% Rich solution after SHEX
h.sol_des_in = (m.sol_poor*h.sol_des_out + m.sol_rich*h.sol_pump_out - m.sol_poor*h.sol_valve_in) / m.sol_rich;
T.sol_des_in = T.sol_pump_out + (h.sol_des_in-h.sol_pump_out)/cp.sol_pump_out;
%% Poor solution after valve Valve
h.sol_abs_in = h.sol_valve_in;
T.sol_abs_in = T.sol_valve_in;
%-------------------------------------------------------------------------%
%% Post processing
% Fluxes over system boundary
Q.cond = m.ref*(h.ref_cond_out-h.ref_cond_in);
Q.evap = m.ref*(h.ref_evap_out-h.ref_evap_in);
Q.abs =  m.sol_rich*h.sol_abs_out - m.ref*h.ref_abs_in - m.sol_poor*h.sol_abs_in;
PP.W_pump = m.sol_rich*(h.sol_pump_out-h.sol_abs_out);
Q.des = m.ref*h.ref_des_out + m.sol_poor*h.sol_des_out - m.sol_rich*h.sol_des_in;
% Heat Exchanger
Q.SHEX = m.sol_poor*(h.sol_des_out-h.sol_valve_in);
Q.RHEX = m.ref*(h.ref_abs_in-h.ref_evap_out);
h.RHEXideal = CoolProp.PropsSI('H','T',T.ref_cond_out,'P',p.evap,'AMMONIA');
eta.RHEX = (h.ref_abs_in-h.ref_evap_out)/(h.RHEXideal-h.ref_evap_out);
eta.SHEX = (h.sol_des_in-h.sol_pump_out)/(h.sol_des_out-h.sol_pump_out);
% COP
PP.COP = Q.evap/(Q.des + PP.W_pump);
% Throttle loss
PP.my_throttle = CoolProp.PropsSI('Q','H',h.ref_valve_in,'P',p.evap,'AMMONIA');
% Circulation
PP.f = m.sol_rich/m.ref;
% Energy balance
PP.energyBalance = Q.des + Q.evap + PP.W_pump + Q.cond + Q.abs;
% Mass balance
PP.massBalance = m.ref + m.sol_poor - m.sol_rich;
%-------------------------------------------------------------------------%
%% Check
% Refrigerant concentrations
if (w.NH3_rich < w.NH3_poor)
    error("w_NH3_rich < w_NH3_poor")
end
if (w.NH3_rich < 0)
    error("w_NH3_rich < 0")
end
if (w.NH3_poor < 0)
    error("w_NH3_poor < 0")
end
if (w.NH3_rich - w.NH3_poor < 0.005)
    error("w_NH3_rich - w_NH3_poor < 0.005")
end
% Mass flows
if (m.ref<0 || m.sol_poor<0 || m.sol_rich<0)
    error("mass flow is negativ")
end
% Energy and mass balance
if (abs(PP.energyBalance) > 1)
    error("Energy is not conserved")
end
if (abs(PP.massBalance) > 1)
    error("Mass is not conserved")
end
%-------------------------------------------------------------------------%
end
