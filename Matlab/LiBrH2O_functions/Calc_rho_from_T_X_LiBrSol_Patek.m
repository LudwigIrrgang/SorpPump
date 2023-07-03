function rho = Calc_rho_from_T_X_LiBrSol_Patek(T,X_mol_LiBr)
% ---------------------------------------------------------------------- %
% Calc_rho_from_T_X_LiBrSol_Patek
% Uses coefficients and formula obtained from Patek 2006
% Input:
%       -   Temperature of Solution T                                   [K]
%       -   Molar Concentration X of LiBr in Solution                   [-]
% Output:
%       -   Molar density rho                                     [mol/m^3]
% ---------------------------------------------------------------------- %
if nargin<2||isempty(X_mol_LiBr),error('Input Argument:Concentration missing');end
if nargin<1||isempty(T),error('Input Argument:Temperature missing');end
%% Constants
T_c = 647.096;              %[K]
rho_c = 17.873;             %[mol/m^3]
RMS = 0.44 * 10^-3;         %[-]
% Table 5
Koef_a = [1.746;4.709];
Koef_t = [0;6];
Koef_m = [1;1];
% Table 12
Koef_beta = [1/3;2/3;5/3;16/3;43/3;110/3];
Koef_alpha = [  1.99274064;
                1.09965342;
                -0.510839303;
                -1.75493479;
                -45.5170352;
                -6.7469445  *   10^5];
%% Calculation
% Calculation of rho_sat
sum = 0;
for i=1:1:6
    sum = sum + Koef_alpha(i)*(1-(T/T_c))^Koef_beta(i);
end
rho_sat = rho_c * (1 + sum);
% Calculation of rho
rho = (1-X_mol_LiBr)*rho_sat+...
    rho_c*(Koef_a(1)*X_mol_LiBr^Koef_m(1)*(T/T_c)^Koef_t(1)+...
    Koef_a(2)*X_mol_LiBr^Koef_m(2)*(T/T_c)^Koef_t(2));
rho = rho*1000;
end
