function s = Calc_s_from_T_X_LiBrSol_Patek(T,X_mol_LiBr)
% ---------------------------------------------------------------------- %
% Calc_s_from_T_X_LiBrSol_Patek
% Uses coefficients and formula obtained from Patek and Klomfar 2006
% DOI: 10.1016/j.ijrefrig.2005.10.007
% Input:
%       -   Temperature of Solution T                                   [K]
%       -   Molar concentration of LiBr in Solution                     [-]
% Output:
%       -   Molar entropy of solution                             [J/mol/K]
% ---------------------------------------------------------------------- %
if nargin<2||isempty(X_mol_LiBr),error('Input Argument:Concentration missing');end
if nargin<1||isempty(T),error('Input Argument:Temperature missing');end
%% Constants
T_c = 647.096;              %[K]
s_c = 79.3933;              %[J/molK]
T_0 = 221;                  %[K]
% Table 8
Koef_a = [  1.53091     *   10^0;
            -4.52564    *   10^0;
            6.98302     *   10^2;
            -2.1666     *   10^4;
            -1.47533    *   10^3;
            8.47012     *   10^-2;
            -6.59523    *   10^0;
            -2.95331    *   10^4;
            9.56314     *   10^-3;
            -1.88679    *   10^-1;
            9.31752     *   10^0;
            5.78104     *   10^0;
            1.38931     *   10^4;
            -1.71762    *   10^4;
            4.15108     *   10^2;
            -5.55647    *   10^4;
            -4.23409    *   10^-3;
            3.05242     *   10^1;
            -1.67620    *   10^0;
            1.48283     *   10^1;
            3.03055     *   10^-3;
            -4.01810    *   10^-2;
            1.49252     *   10^-1;
            2.59240     *   10^0;
            -1.77421    *   10^-1;
            -6.99650    *   10^-5;
            6.05007     *   10^-4;
            -1.65228    *   10^-3;
            1.22966     *   10^-3];
Koef_t = [0;0;0;0;0;1;1;1;2;2;2;2;2;2;2;2;3;3;3;3;4;4;4;4;4;5;5;5;5];
Koef_n = [0;1;6;6;2;0;0;4;0;0;4;0;4;5;2;5;0;4;0;1;0;2;4;7;1;0;1;2;3];
Koef_m = [1;1;2;3;6;1;3;5;1;2;2;4;5;5;6;6;1;3;5;7;1;1;1;2;3;1;1;1;1];
% Table 15
Koef_beta = [1/3;1;8/3;8];
Koef_alpha = [  -3.34112    *   10^-1;
                -8.47987    *   10^-1;
                -9.11980    *   10^-1;
                -1.64046    *   10^0];
%% Calculation
% Calculation of s_sat
sum = 0;
for i=1:1:4
    sum = sum + Koef_alpha(i)*(1-(T/T_c))^Koef_beta(i);
end
s_sat = s_c * (1 + sum);
% Calculation of s
factors = zeros(29,1);
a = 0;
b = 0;
c = 0;
d = 0;
e = 0;
f = 0;
for i=1:1:29
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
s = (1-X_mol_LiBr)*s_sat + s_c*(a*(T_c/(T-T_0))^5 + b*(T_c/(T-T_0))^4 + c*(T_c/(T-T_0))^3 + d*(T_c/(T-T_0))^2 + e*(T_c/(T-T_0))^1 + f);
end
