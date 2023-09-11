%% Script for Use of AHP Model 
% Defines the necessary input variables for the AHP model and calls it
% Use this script to investigate certain operational points
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
The distibution of this script and all incoperated models is not permitted
without permission of the owner.
%}
% ----------------------------------------------------------------------- %
%% CoolProp, Refrop, Functions
% Coolprop (update path!)
addpath('M:\Benutzer\Irrgang\03_Numerische Untersuchungen\Numerical Model AHP\CoolProp')
% Refprop (upadate path if necessary)
addpath('C:\Program Files (x86)\REFPROP');
% Functions
addpath('.\H2ONH3_functions')
addpath('.\LiBrH2O_functions')
% Internal Cycle Models
addpath('.\Internal Cycle Models')
% Cycle plots
addpath('.\Cycle Plot')
% ----------------------------------------------------------------------- %
%% Boundary conditions/ Input variables
% ----------------------------------------------------------------------- %
% Efficiencies
    % Pump efficiency
    eta.pump = 1;
% Temperatures
    % Heat source temperature
    T.ext_des_in = 160+273.15;
    % Production temperature
    T.ext_evap_out = 8+273.15;
    T.ext_evap_in = 15+273.15;
    % Heat sink temperature (Parallel cooling - same temperature at both HXs)
    T.ext_abs_in = 30+273.15;
    T.ext_cond_in = 30+273.15;
% Model setup
    % Use 'Q_evap' for cold output/ 'Q_des' for heat input
    s.requirement = 'Q_evap';
    Q.dec = 10*10^3;
    % Configuration
    % Use 'Base', 'Double_effect', 'Double_lift'
    s.configuration = 'Base';
    % Refrigerant 
    % Use 'AMMONIA' or 'WATER'
    s.refrigerant = 'AMMONIA';
    % Cycle plot
    % Use 1 if you want to plot the pT diagram
    s.plot_Duhring = 1;
    % Minimal difference in concentration between solution flows
    s.delta_w = 0.05;
% Approach temperatures
    HX.T_PP_evap =      5;       % Evaporation approach temperature
    HX.T_PP_abs =       5;       % Absorber approach temperature
    HX.T_PP_des =       5;       % Desorber approach temperature
    HX.T_PP_cond =      5;       % Condenser approach temperature
    HX.T_PP_SHEX =      5;       % Solution heat exchanger approach temperature
    HX.T_PP_SHEXI =     5;       % Solution heat exchanger approach temperature
    HX.T_PP_RHEX =      5;       % Refrigerant heat exchanger approach temperature
    HX.T_PP_cond_int =  5;       % Internal condenser approach temperature
    HX.SC_cond =        5;       % Subcooling at condenser
    HX.dT_ref_des =     5;       % Difference between solution and refrigerant desober exit temperatures
    HX.dT_ref_desI =    5;       % Difference between solution and refrigerant desoberI exit temperatures
% ----------------------------------------------------------------------- %
%% Run AHP model
[T, p, h, m, w, eta, Q, PP, s] = AHP(T,HX,Q,s,eta)
