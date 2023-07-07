function wNH3 = NH3inVapor_Calc_X_PT(Pressure,Temp)
% ---------------------------------------------------------------------- %
% NH3inVapor_Calc_X_PT
% This function calculates the concentration of the saturated NH3-H2O
% vapor from pressure and temperature
% Make sure that REFPROP can be called correctly
% Input units: Pressure [Pa], Tempreature [K]
% ---------------------------------------------------------------------- %
if nargin<2||isempty(Pressure),error('Input Argument:Pressure missing');end
if nargin<1||isempty(Temp),error('Input Argument:Temperature missing');end
%% Calculation
i = 0;
a = 1;
b = 0;
try
    maxPressure = refpropm('P','T',Temp,'Q',1,'AMMONIA','WATER',[a b]);
catch
    warning('Refprop failed in two phase region - only ammonia in vapor')
    wNH3 = 1; % only ammonia present in vapor
    return
end
if(maxPressure<Pressure)
    warning('Pressure not in two phase region')
    wNH3 = 1; % only water present in vapor
else
    f=@(x)Pressure-refpropm('P','T',Temp,'Q',1,'AMMONIA','WATER',[x 1-x]);
    while abs(f(a))>1 % corresponds to 0.1 kPa
        p = refpropm('P','T',Temp,'Q',1,'AMMONIA','WATER',[a 1-a]);
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
    wNH3 = a;
end
