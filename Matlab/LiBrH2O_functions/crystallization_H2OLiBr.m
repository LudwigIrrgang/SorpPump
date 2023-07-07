function [CritValue] = crystallization_H2OLiBr(decider,value)
% ---------------------------------------------------------------------- %
% crystallization
% Uses data and formulations from Boryta 1970
% Calculates critical Temperatures or concentrations for given value
% Use "T" if Input is temperature as decider value
% Use "w" if Input is mass fraction of LiBr in solution as decider value
% ---------------------------------------------------------------------- %
%% Computation
T_cr = 0;
w_cr = 0;
if decider == "w"
    % Parameters (from Albers report 2019, Boryta 1970)
    parameter_T = [ 42.90198341384762;
                    34.67510890651030;
                    31.30778644395644;
                    2.99859601946791;
                    -19.36781324384540;
                    -4.88856108511827;
                    4.61433775768846;
                    1.80636830673333    ];
    % Define working variable
    w = value;
    w_r = (w-0.64794)/0.044858;
    % Calculate critical values
    for i=1:1:8
        T_cr = T_cr + parameter_T(i)*(w_r^(i-1));
    end
    CritValue = T_cr;
elseif decider == "T"
    % Parameters (from Albers report 2019, Boryta 1970)
    parameter_w = [ 0.66136507494441;
                    0.02262634534253;
                    -0.02216522722755;
                    0.05134156572205;
                    0.00034455919818;
                    -0.03628931060739;
                    0.00252166562759;
                    0.00796985214167    ];
    % Define working variable
    T = value - 273.15;
    T_r = (T-54.793)/33.111;
    % Calculate critical values
    for i=1:1:8
        w_cr = w_cr + parameter_w(i)*(T_r^(i-1));
    end
    CritValue = w_cr;
else
    error("Crystallization function failed - Decider value not defined - use T or w")
end
end
