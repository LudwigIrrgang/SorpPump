function [T_abs_in] = NH3inSolution_Calc_state_valve_exit(T_sat, p_evap, h_abs_in, m_poor, w_NH3_poor, stepsize)
% ---------------------------------------------------------------------- %
% Calc_state_valve_exit
% Since the solution is depressurized beyond saturation some refrigerant
% has to be evaporated
% The temperature is decreases until mass and energy equilibrium are
% reached
% ---------------------------------------------------------------------- %
if nargin<6||isempty(T_sat),error('Input Argument:T_sat missing');end
if nargin<5||isempty(p_evap),error('Input Argument:p_evap missing');end
if nargin<4||isempty(h_abs_in),error('Input Argument:h_abs_in missing');end
if nargin<3||isempty(m_poor),error('Input Argument:m_poor missing');end
if nargin<2||isempty(w_NH3_poor),error('Input Argument:w_NH3_poor missing');end
if nargin<1||isempty(stepsize),error('Input Argument:stepsize missing');end
%% Calculation
% saturation as starting point
w_NH3 = NH3inSolution_Calc_X_PT(p_evap/1000,T_sat);
h_mixture = refpropm('H','T',T_sat,'Q',0,'AMMONIA','WATER',[w_NH3 (1-w_NH3)]);
T_abs_in = T_sat - stepsize;
% loop
while(h_mixture>h_abs_in) % Decrease the temerature until energy equilibrium is achieved
    w_NH3 = NH3inSolution_Calc_X_PT(p_evap/1000,T_des_in);
    w_NH3_vapor = NH3inVapor_Calc_X_PT(p_cond/1000,T_des_in);
    h_sol = refpropm('H','T',T_abs_in,'Q',0,'AMMONIA','WATER',[w_NH3 (1-w_NH3)]);
    h_vapor = refpropm('H','T',T_abs_in,'Q',1,'AMMONIA','WATER',[w_NH3_vapor (1-w_NH3_vapor)]);
    A = [   w_NH3,      w_NH3_vapor;
            1,          1           ];
    b = [   m_poor*w_NH3_rich;     m_poor];
    x = A\b;
    m_sol = x(1);
    m_vapor = x(2);
    h_mixture = (m_sol*h_sol + m_vapor*h_vapor)/(m_sol+m_vapor);
    if(h_mixture>h_abs_in)
        T_abs_in = T_abs_in - stepsize;
    end
end
end