function [T_abs_in] = Calc_state_valve_exit(T_sat, p_evap, h_abs_in, m_poor, w_H2O_poor, stepsize)
% ---------------------------------------------------------------------- %
% Calc_state_valve_exit
% Input:
%       -   Saturation temperature of Solution T                    [K]
%       -   Evaporation pressure                                    [Pa]
%       -   Enthalpy of solution absorber inlet                     [J/kg]
%       -   Poor solution mass flow rate                            [kg/s]
%       -   Poor solution H2O concentration                         [-]
%       -   Stepzise                                                [-]
% Output:
%       -   Solution temperature absorber inlet                     [K]
% ---------------------------------------------------------------------- %
% Since the solution is depressurized beyond saturation some refrigerant
% has to be evaporated
% The solution temperature is decreases until mass and energy equilibrium
% are reached
% ---------------------------------------------------------------------- %
if nargin<6||isempty(T_sat),error('Input Argument:T_sat missing');end
if nargin<5||isempty(p_evap),error('Input Argument:p_evap missing');end
if nargin<4||isempty(h_abs_in),error('Input Argument:h_abs_in missing');end
if nargin<3||isempty(m_poor),error('Input Argument:m_poor missing');end
if nargin<2||isempty(w_H2O_poor),error('Input Argument:w_H2O_poor missing');end
if nargin<1||isempty(stepsize),error('Input Argument:stepsize missing');end
%% Definition of constants
% ------------------------Necessary Constants---------------------------- %
M_LiBr = 0.08685;               %[kg/mol]
M_H2O = 0.018015268;            %[kg/mol]
%% Calculation
% saturation as starting point
x_LiBr = Calc_X_from_T_p_satLiBrSol_Patek(T_sat,p_evap);
h_mixture = Calc_h_from_T_X_LiBrSol_Patek(T_sat,x_LiBr);
T_abs_in = T_sat - stepsize;
% loop
while(h_mixture>h_abs_in) % Decrease the temperature until energy equilibrium is achieved
    x_LiBr = Calc_X_from_T_p_satLiBrSol_Patek(T_abs_in,p_evap);
    w_LiBr = x_LiBr*M_LiBr/ (x_LiBr*M_LiBr + (1-x_LiBr)*M_H2O);
    w_H2O = 1 - w_LiBr;
    m_H2O = (w_H2O_poor - w_H2O)*m_poor;
    m_sol = m_poor - m_H2O;
    h_sol_mol = Calc_h_from_T_X_LiBrSol_Patek(T_abs_in,x_LiBr);
    h_sol = h_sol_mol/(x_LiBr*M_LiBr + (1-x_LiBr)*M_H2O);
    h_H2O = CoolProp.PropsSI('H','T',T_abs_in,'P',p_evap,'Water');
    h_mixture = (m_sol*h_sol + m_H2O*h_H2O)/(m_sol+m_H2O);
    if(h_mixture<h_des_in)
        T_abs_in = T_abs_in - stepsize;
    end
end
end
