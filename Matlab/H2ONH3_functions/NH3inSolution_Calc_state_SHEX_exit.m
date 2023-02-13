function [T_des_in] = NH3inSolution_Calc_state_SHEX_exit(T_sat, p_cond, h_des_in, m_rich, w_NH3_rich, stepsize)
% ---------------------------------------------------------------------- %
% Calc_state_SHEX_exit
% Since the solution is heated above saturation, some refrigerant must be
% evaporated
% The solution temperature is increased until mass and energy equilibrium
% are reached
% ---------------------------------------------------------------------- %
if nargin<6||isempty(T_sat),error('Input Argument:T_sat missing');end
if nargin<5||isempty(p_cond),error('Input Argument:p_cond missing');end
if nargin<4||isempty(h_des_in),error('Input Argument:h_des_in missing');end
if nargin<3||isempty(m_rich),error('Input Argument:m_rich missing');end
if nargin<2||isempty(w_NH3_rich),error('Input Argument:w_NH3_rich missing');end
if nargin<1||isempty(stepsize),error('Input Argument:stepsize missing');end
%% Calculation
% saturation as starting point
w_NH3 = NH3inSolution_Calc_X_PT(p_cond/1000,T_sat);
h_mixture = refpropm('H','T',T_sat,'Q',0,'AMMONIA','WATER',[w_NH3 (1-w_NH3)]); % enthalpy of satureated solution 
T_des_in = T_sat;
% loop
while(h_mixture<h_des_in) % Increase the temerature until energy equilibrium is achieved
    w_NH3 = NH3inSolution_Calc_X_PT(p_cond/1000,T_des_in);
    w_NH3_vapor = NH3inVapor_Calc_X_PT(p_cond/1000,T_des_in);
    %w_NH3_vapor = 1;
    h_sol = refpropm('H','T',T_des_in,'Q',0,'AMMONIA','WATER',[w_NH3 (1-w_NH3)]);
    h_sol_p = refpropm('H','T',T_des_in,'P',p_cond/1000,'AMMONIA','WATER',[w_NH3 (1-w_NH3)]);
    h_vapor = refpropm('H','T',T_des_in,'Q',1,'AMMONIA','WATER',[w_NH3_vapor (1-w_NH3_vapor)]);
    h_vapor_p = refpropm('H','T',T_des_in,'P',p_cond/1000,'AMMONIA','WATER',[w_NH3_vapor (1-w_NH3_vapor)]);
    A = [   w_NH3,      w_NH3_vapor;
            1,          1           ];
    b = [   m_rich*w_NH3_rich;     m_rich];
    x = A\b;
    m_sol = x(1);
    m_vapor = x(2);
    h_mixture = (m_sol*h_sol + m_vapor*h_vapor)/(m_sol+m_vapor);
    if(h_mixture<h_des_in)
        T_des_in = T_des_in + stepsize;
    end
end
end