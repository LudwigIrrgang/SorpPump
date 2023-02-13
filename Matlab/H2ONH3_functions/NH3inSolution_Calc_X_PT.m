function wNH3 = NH3inSolution_Calc_X_PT(Pressure,Temp)
% ---------------------------------------------------------------------- %
% NH3inSolution_Calc_X_PT
% This function calculates the concentration of the saturated NH3-H2O
% solution from pressure and temperature
% Input units: Pressure [Pa], Tempreature [K]
% ---------------------------------------------------------------------- %
if nargin<2||isempty(Pressure),error('Input Argument:Pressure missing');end
if nargin<1||isempty(Temp),error('Input Argument:Temperature missing');end
%% Calculation
i = 0;
a = 0.999;
b = 0.001;
try
    minPressure = refpropm('P','T',Temp,'Q',0,'WATER','AMMONIA',[a b]);
catch
    warning('Refprop failed in two phase region - only water in solution')
    wNH3 = 0; % only water present in solution
    return
end
if(minPressure>Pressure)
    warning('Pressure not in two phase region')
    wNH3 = 0; % only water present in solution
else
    f=@(x)Pressure-refpropm('P','T',Temp,'Q',0,'WATER','AMMONIA',[1-x x]);
    while abs(f(b))>0.1 % corresponds to 0.1 kPa
        if f((a+b)/2)<0 % half splitting technique
            a=(a+b)/2;
        else
            b=(a+b)/2;
        end
        i = i+1;
        if (i>1000)
            error('Concentration can not be calculated - loop exceeded limits')
        end
    end
    wNH3 = b;
end
