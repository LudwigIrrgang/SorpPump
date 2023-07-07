function [T, p, h, m, w, eta, Q, PP, s] = AHP(T,HX,Q,s,eta)
%% Function Absorption Heat Pump
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
% Calculates internal and external state points of:
% - Single effect
% - Double effect (parallel)
% - Double lift
% Absorption chiller with the following working fluids:
% - Lithium bromide/ water
% - Water/ ammonia
% Make sure that all necessary functions can be called correctly
% See run_AHP.m for example code
% ----------------------------------------------------------------------- %
%% Initialize Structs
h(1) = struct();
m(1) = struct();
p(1) = struct();
% ----------------------------------------------------------------------- %
%% Define internal model input 
% Evaporator
T.evap = T.ext_evap_out - HX.T_PP_evap;
% Absorber
T.sol_abs_out = T.ext_abs_in + HX.T_PP_abs;
% Desorber
T.sol_des_out = T.ext_des_in - HX.T_PP_des;
% Condenser
T.cond = T.ext_cond_in + HX.T_PP_cond + HX.SC_cond;
% Specific values for configuration
switch s.configuration
    case 'Double_effect'
        % Internal condenser
        T.cond_int = T.cond + HX.dT_ref_des + HX.T_PP_cond_int;
        % Upper pressure desorber
        T.sol_des_outI = T.sol_des_out;
        % Middle pressur desorber
        T.sol_des_out = T.cond_int - HX.T_PP_cond_int;
    case 'Double_lift'
        % Upper pressure desorber
        T.sol_des_outI = T.sol_des_out;
        % Middle pressure absorber
        T.sol_abs_outI = T.sol_abs_out;
    otherwise
end
% ----------------------------------------------------------------------- %
%% Calculation
switch s.configuration
% ----------------------------------------------------------------------- %
% -------------------------- Base Model --------------------------------- %
% ----------------------------------------------------------------------- %
    case 'Base'
        %% Calculate absorption machine internals
        switch s.refrigerant
            case 'WATER'
                [T, p, h,  m, w, eta, Q, PP, s] = base_model_H2OLiBr(T,p,h,m,eta,Q,HX,s);
            case 'AMMONIA'
                [T, p, h,  m, w, eta, Q, PP, s] = base_model_NH3H2O(T,p,h,m,eta,Q,HX,s);
            otherwise
                error('Refrigerant not available, use WATER or AMMONIA')
        end
% ----------------------------------------------------------------------- %
% ---------------------- Double Effect Model ---------------------------- %
% ----------------------------------------------------------------------- %
    case 'Double_effect'
        %% Calculate absorption machine internals
        switch s.refrigerant
            case 'WATER'
                [T, p, h, m, w, eta, Q, PP, s] = doubleEffect_model_H2OLiBr(T,p,h,m,eta,Q,HX,s);
            case 'AMMONIA'
                [T, p, h, m, w, eta, Q, PP, s] = doubleEffect_model_NH3H2O(T,p,h,m,eta,Q,HX,s);
            otherwise
                error('Refrigerant not available, use WATER or AMMONIA')
        end
% % ----------------------------------------------------------------------- %
% % ------------------------ Double Lift Model ---------------------------- %
% % ----------------------------------------------------------------------- %
    case 'Double_lift'
        %% Calculate absorption machine internals
        switch s.refrigerant
            case 'WATER'
                [T, p, h, m, w, eta, Q, PP, s] = doubleLift_model_H2OLiBr(T,p,h,m,eta,Q,HX,s);
            case 'AMMONIA'
                [T, p, h, m, w, eta, Q, PP, s] = doubleLift_model_NH3H2O(T,p,h,m,eta,Q,HX,s);
            otherwise
                error('Refrigerant not available, use WATER or AMMONIA')
        end
    otherwise 
        error('Configuration not available')
end
%% Duhring Plot
if s.plot_Duhring == 1 && s.refrigerant == "WATER"
    LiBr_pTDiagram(T, p, s);
elseif s.plot_Duhring == 1 && s.refrigerant == "AMMONIA"
    NH3_pTDiagram(T, p, s);
end
end
