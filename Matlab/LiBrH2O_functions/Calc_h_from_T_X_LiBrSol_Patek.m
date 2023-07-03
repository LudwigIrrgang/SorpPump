function h = Calc_h_from_T_X_LiBrSol_Patek(T,X_mol_LiBr)
% ---------------------------------------------------------------------- %
% Calc_h_from_T_X_LiBrSol_Patek
% Uses coefficients and formula obtained from Patek 2006
% Input:
%       -   Temperature of Solution T                                   [K]
%       -   Molar Concentration X of LiBr in Solution                   [-]
% Output:
%       -   Molar enthalpy h of solution                            [J/mol]
% ---------------------------------------------------------------------- %
if nargin<2||isempty(X_mol_LiBr),error('Input Argument:Concentration missing');end
if nargin<1||isempty(T),error('Input Argument:Temperature missing');end
%% Constants
T_c = 647.096;              %[K]
h_c = 37548.5;              %[J/mol]
T_0 = 221;                  %[K]
% Table 7
Koef_a = [  2.27431     *   10^0;
            -7.99511    *   10^0;
            3.85239     *   10^2;
            -1.63940    *   10^4;
            -4.22562    *   10^2;          
            1.13314     *   10^-1;
            -8.33474    *   10^0;
            -1.73833    *   10^4;
            6.49763     *   10^0;
            3.24552     *   10^3;
            -1.34643    *   10^4;
            3.99322     *   10^4;
            -2.58877    *   10^5;
            -1.93046    *   10^-3;
            2.80616     *   10^0;
            -4.04479    *   10^1;
            1.45342     *   10^2;
            -2.74873    *   10^0;
            -4.49743    *   10^2;
            -1.21794    *   10^1;
            -5.83739    *   10^-3;
            2.33910     *   10^-1;
            3.41888     *   10^-1;
            8.85259     *   10^0;
            -1.78731    *   10^1;
            7.35179     *   10^-2;
            -1.79430    *   10^-4;
            1.84261     *   10^-3;
            -6.24282    *   10^-3;
            6.84765     *   10^-3];
Koef_t = [0;0;0;0;0;1;1;1;2;2;2;2;2;3;3;3;3;3;3;3;4;4;4;4;4;4;5;5;5;5];
Koef_n = [0;1;6;6;2;0;0;4;0;4;5;5;6;0;3;5;7;0;3;1;0;4;2;6;7;0;0;1;2;3];
Koef_m = [1;1;2;3;6;1;3;5;4;5;5;6;6;1;2;2;2;5;6;7;1;1;2;2;2;3;1;1;1;1];
% Table 14
Koef_beta = [1/3;2/3;5/6;21/6];
Koef_alpha = [  -4.37196    *   10^-1;
                3.03440     *   10^-1;
                -1.29582    *   10^0;
                -1.76410    *   10^-1];
%% Calculation
% Calculation of h_sat
sum = 0;
for i=1:1:4
    sum = sum + Koef_alpha(i)*(1-(T/T_c))^Koef_beta(i);
end
h_sat = h_c * (1 + sum);
% Calculation of h
factors = zeros(30,1);
a = 0;
b = 0;
c = 0;
d = 0;
e = 0;
f = 0;
for i=1:1:30
   factors(i) =  Koef_a(i)*X_mol_LiBr^Koef_m(i)*(0.4-X_mol_LiBr)^Koef_n(i);
   if Koef_t(i) == 0
       f = f + factors(i);
   elseif Koef_t(i) == 1
       e = e + factors(i);
   elseif Koef_t(i) == 2
       d = d + factors(i);
   elseif Koef_t(i) == 3
       c = c + factors(i);
   elseif Koef_t(i) == 4
       b = b + factors(i);
   elseif Koef_t(i) == 5
       a = a + factors(i);
   end
end
h = (1-X_mol_LiBr)*h_sat + h_c*(a*(T_c/(T-T_0))^5 + b*(T_c/(T-T_0))^4 + c*(T_c/(T-T_0))^3 + d*(T_c/(T-T_0))^2 + e*(T_c/(T-T_0))^1 + f);
end
