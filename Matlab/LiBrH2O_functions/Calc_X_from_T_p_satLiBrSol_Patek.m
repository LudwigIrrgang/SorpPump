function X = Calc_X_from_T_p_satLiBrSol_Patek(T,p)
% ---------------------------------------------------------------------- %
% Calc_X_from_T_p_LiBrSol_Patek
% Uses coefficients and formula obtained from Patek and Klomfar 2006
% DOI: 10.1016/j.ijrefrig.2005.10.007
% Input:
%       -   Temperature of LiBR solution                                [K]
%       -   Pressure of saturated solution                             [Pa]
% Output:
%       -   Molar concentration of LiBr in saturated solution           [-]
% ---------------------------------------------------------------------- %
if nargin<2||isempty(p),error('Input Argument:Pressure missing');end
if nargin<1||isempty(T),error('Input Argument:Temperature missing');end
%% Constants
T_c = 647.096;              %[K]
p_c = 22.064 * 10^6;        %[Pa]
% Table 4
Koef_a = [  -2.41303    *   10^2;
            1.91750     *   10^7;
            -1.75521    *   10^8;
            3.25430     *   10^7;
            3.92571     *   10^2;
            -2.12626    *   10^3;
            1.85127     *   10^8;
            1.91216     *   10^3];
Koef_t = [0;0;0;0;1;1;1;1];
Koef_n = [0;5;6;3;0;2;6;0];
Koef_m = [3;4;4;8;1;1;4;6];
% Table 11
Koef_beta = [1.0;1.5;3.0;3.5;4.0;7.5];
Koef_alpha = [  -7.85951783;
                1.84408259;
                -11.7866497;
                22.6807411;
                -15.9618719;
                1.80122502];
%% Calculation
% Calculate T_sat
funT = @(T)p_c*exp((T_c/T)*...
    (Koef_alpha(1)*(1-(T/T_c))^Koef_beta(1)+...
    Koef_alpha(2)*(1-(T/T_c))^Koef_beta(2)+...
    Koef_alpha(3)*(1-(T/T_c))^Koef_beta(3)+...
    Koef_alpha(4)*(1-(T/T_c))^Koef_beta(4)+...
    Koef_alpha(5)*(1-(T/T_c))^Koef_beta(5)+...
    Koef_alpha(6)*(1-(T/T_c))^Koef_beta(6)))-p;
T_sat = fzero(funT,400);

T_R = T/T_c;
funX = @(X)T-(Koef_a(1)*X^Koef_m(1)*(0.4-X)^Koef_n(1)*(T_R)^Koef_t(1)...
    +Koef_a(2)*X^Koef_m(2)*(0.4-X)^Koef_n(2)*(T_R)^Koef_t(2)...
    +Koef_a(3)*X^Koef_m(3)*(0.4-X)^Koef_n(3)*(T_R)^Koef_t(3)...
    +Koef_a(4)*X^Koef_m(4)*(0.4-X)^Koef_n(4)*(T_R)^Koef_t(4)...
    +Koef_a(5)*X^Koef_m(5)*(0.4-X)^Koef_n(5)*(T_R)^Koef_t(5)...
    +Koef_a(6)*X^Koef_m(6)*(0.4-X)^Koef_n(6)*(T_R)^Koef_t(6)...
    +Koef_a(7)*X^Koef_m(7)*(0.4-X)^Koef_n(7)*(T_R)^Koef_t(7)...
    +Koef_a(8)*X^Koef_m(8)*(0.4-X)^Koef_n(8)*(T_R)^Koef_t(8))-T_sat;
X = fzero(funX,0);
end
