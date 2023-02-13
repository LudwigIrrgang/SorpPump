function [T_des_in] = Calc_state_SHEX_exit(T_sat, p_cond, h_des_in, m_rich, w_H2O_rich, stepsize)
% ---------------------------------------------------------------------- %
% Calc_state_SHEX_exit
% ---------------------------------------------------------------------- %
%{
Author  : Ludwig Irrgang
Date    : 25.06.2022
Copyright information:
Ludwig Irrgang
Lehrstuhl für Energiesysteme
TUM School of Engineering and Design
Technische Universität München
Boltzmannstr. 15 
85748 Garching b. München
ludwig.irrgang@tum.de
%}
% ---------------------------------------------------------------------- %
% Input:
%       -   Saturation temperature of Solution T                    [K]
%       -   Condensation pressure                                   [Pa]
%       -   Enthalpy of solution decorber inlet                     [J/kg]
%       -   Rich solution mass flow rate                            [kg/s]
%       -   Rich solution H2O concentration                         [-]
%       -   Stepzise                                                [-]
% Output:
%       -   Solution temperature desorber inlet                     [K]
% ---------------------------------------------------------------------- %
% Since the solution is heated above saturation, some refrigerant must be
% evaporated
% The solution temperature is increased until mass and energy equilibrium
% are reached
% ---------------------------------------------------------------------- %
if nargin<6||isempty(T_sat),error('Input Argument:T_sat missing');end
if nargin<5||isempty(p_cond),error('Input Argument:p_cond missing');end
if nargin<4||isempty(h_des_in),error('Input Argument:h_des_in missing');end
if nargin<3||isempty(m_rich),error('Input Argument:w_H2O_rich missing');end
if nargin<2||isempty(w_H2O_rich),error('Input Argument:w_H2O_rich missing');end
if nargin<1||isempty(stepsize),error('Input Argument:stepsize missing');end
%% Definition of constants
% ------------------------Necessary Constants---------------------------- %
M_LiBr = 0.08685;               %[kg/mol]
M_H2O = 0.018015268;            %[kg/mol]
%% Calculation
% saturation as starting point
x_LiBr = Calc_X_from_T_p_satLiBrSol_Patek(T_sat,p_cond);
h_mixture = Calc_h_from_T_X_LiBrSol_Patek(T_sat,x_LiBr);
T_des_in = T_sat + stepsize;
% loop
while(h_mixture<h_des_in) % Increase the temerature until energy equilibrium is achieved
    x_LiBr = Calc_X_from_T_p_satLiBrSol_Patek(T_des_in,p_cond);
    w_LiBr = x_LiBr*M_LiBr/ (x_LiBr*M_LiBr + (1-x_LiBr)*M_H2O);
    w_H2O = 1 - w_LiBr;
    m_H2O = (w_H2O_rich - w_H2O)*m_rich;
    m_sol = m_rich - m_H2O;
    h_sol_mol = Calc_h_from_T_X_LiBrSol_Patek(T_des_in,x_LiBr);
    h_sol = h_sol_mol/(x_LiBr*M_LiBr + (1-x_LiBr)*M_H2O);
    h_H2O = CoolProp.PropsSI('H','T',T_des_in,'P',p_cond,'Water');
    h_mixture = (m_sol*h_sol + m_H2O*h_H2O)/(m_sol+m_H2O);
    if(h_mixture<h_des_in)
        T_des_in = T_des_in + stepsize;
    end
end
end