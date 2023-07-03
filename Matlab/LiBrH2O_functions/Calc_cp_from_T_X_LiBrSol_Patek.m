function cp = Calc_cp_from_T_X_LiBrSol_Patek(T,X_mol_LiBr)
% ---------------------------------------------------------------------- %
% Calc_cp_from_T_X_LiBrSol_Patek
% Uses coefficients and formula obtained from Patek 2006
% Input:
%       -   Temperature of Solution T                                   [K]
%       -   Molar Concentration X of LiBr in Solution                   [-]
% Output:
%       -   Molar heat capacity                                   [J/mol/K]
% ---------------------------------------------------------------------- %
if nargin<2||isempty(X_mol_LiBr),error('Input Argument:Concentration missing');end
if nargin<1||isempty(T),error('Input Argument:Temperature missing');end
%% Constants
cp_t = 76.0226;             %[J/mol/K]
T_c = 647.096;              %[K]
T_0 = 221;                  %[K]
RMS = 0.98 * 10^-3;         %[-]
T_t = 273.16;               %[K]
% Table 6
Koef_a = [-1.42094  *   10^1;
          4.04943   *   10^1;
          1.11135   *   10^2;
          2.29980   *   10^2;
          1.34526   *   10^3;
          -1.41010  *   10^-2;
          1.24977   *   10^-2;
          -6.83209  *   10^-4];
Koef_t = [0;0;0;0;0;2;3;4];
Koef_n = [0;0;1;2;3;0;3;2];
Koef_m = [2;3;3;3;3;2;1;1];
% Table 13
Koef_beta = [0;2;3;6;34];
Koef_gamma = [0;2;3;5;0];
Koef_alpha = [1.38801;
              -2.95318
              3.18721;
              -0.645473;
              9.18946 * 10^5];
%% Calculation
% Calculation of cp_sat
sum = 0;
for i=1:1:5
    sum = sum + Koef_alpha(i)*(1-(T/T_c))^Koef_beta(i)*(T/T_t)^Koef_gamma(i);
end
cp_sat = cp_t*sum;
% Calculation of cp
cp = (1-X_mol_LiBr)*cp_sat+...
    cp_t*(Koef_a(1)*X_mol_LiBr^Koef_m(1)*(0.4 - X_mol_LiBr)^Koef_n(1)*(T_c/(T-T_0))^Koef_t(1)+...
    Koef_a(2)*X_mol_LiBr^Koef_m(2)*(0.4 - X_mol_LiBr)^Koef_n(2)*(T_c/(T-T_0))^Koef_t(2)+...
    Koef_a(3)*X_mol_LiBr^Koef_m(3)*(0.4 - X_mol_LiBr)^Koef_n(3)*(T_c/(T-T_0))^Koef_t(3)+...
    Koef_a(4)*X_mol_LiBr^Koef_m(4)*(0.4 - X_mol_LiBr)^Koef_n(4)*(T_c/(T-T_0))^Koef_t(4)+...
    Koef_a(5)*X_mol_LiBr^Koef_m(5)*(0.4 - X_mol_LiBr)^Koef_n(5)*(T_c/(T-T_0))^Koef_t(5)+...
    Koef_a(6)*X_mol_LiBr^Koef_m(6)*(0.4 - X_mol_LiBr)^Koef_n(6)*(T_c/(T-T_0))^Koef_t(6)+...
    Koef_a(7)*X_mol_LiBr^Koef_m(7)*(0.4 - X_mol_LiBr)^Koef_n(7)*(T_c/(T-T_0))^Koef_t(7)+...
    Koef_a(8)*X_mol_LiBr^Koef_m(8)*(0.4 - X_mol_LiBr)^Koef_n(8)*(T_c/(T-T_0))^Koef_t(8));
end
