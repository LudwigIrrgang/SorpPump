function p_vap = Calc_p_from_T_X_LiBrSol_Patek(T,X_mol_LiBr)
% ---------------------------------------------------------------------- %
% Calc_p_from_T_X_LiBrSol_Patek
% Uses coefficients and formula obtained from Patek 2006
% Input:
%       -   Temperature of Solution T                                   [K]
%       -   Molar Concentration X of LiBr in Solution                   [-]
% Output:
%       -   Pressure of solution at saturation                         [Pa]
% ---------------------------------------------------------------------- %
if nargin<2||isempty(X_mol_LiBr),error('Input Argument:Concentration missing');end
if nargin<1||isempty(T),error('Input Argument:Temperature missing');end
%% Constants
T_c = 647.096;              %[K]
p_c = 22.064e6;             %[Pa]
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
% Calculation of Teta
sum = 0;
for i=1:1:8
   sum = sum +  Koef_a(i)*X_mol_LiBr^Koef_m(i)*(0.4-X_mol_LiBr)^Koef_n(i)*(T/T_c)^Koef_t(i);
end
Teta = T - sum;
% Calculation of vapor pressure p_vap
sum = 0;
for i=1:1:6
    sum = sum + Koef_alpha(i)*(1-(Teta/T_c))^Koef_beta(i);
end
p_vap = p_c*exp((T_c/Teta)*sum);
end
