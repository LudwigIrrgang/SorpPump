function checkForViolation_H2OLiBr(w,T,position)
% ---------------------------------------------------------------------- %
% function checkForViolation
% ---------------------------------------------------------------------- %
%{
Author  : Ludwig Irrgang
Date    : 25.06.2022
Copyright information:
Ludwig Irrgang
Lehrstuhl für Energiesysteme
TUM School of Engineering and Design
Technische Universität München
Boltzmannstr. 15 
85748 Garching b. München
ludwig.irrgang@tum.de
%}
% ---------------------------------------------------------------------- %
%% Description
% Checks for violation of accepted values for solution state functions
% Gives no output in case of no violation
% Throws error if boundaries violated
% Input:    w : mass fraction of LiBr in solution as decider value
%           T : temperature of solution
%% Computation
w_cr = crystallization_H2OLiBr("T",T);
T_cr = crystallization_H2OLiBr("w",w);
% Check for Patek boundaries (273K - 500K)
if T>500
    error("Temperature of LiBr solution ("+T+"K) above Patek max. tempreture of 500 K at point: "+position)
elseif T<=273.15
    error("Temperature of LiBr solution ("+T+"K) below Patek min. tempreture of 273.15 K at point: "+position)
end
% Check for Boryta boundaries (273K - 374K)
if T>374
    %warning("Crystallization can not be checked because Temperature ("+T+"K) of LiBr solution above Boryta max. tempreture of 374 K at point: "+position)
elseif T<=273.15
    %warning("Crystallization can not be checked because Temperature of LiBr solution ("+T+"K) below min. tempreture of 273.15 K at point: "+position)
end
if T<=374 && T>=273.15 && w>0.57
    if w>w_cr
        error("Crystallization Error! Mass fraction of LiBr in solution("+w+") above max. value of "+w_cr+" at position: "+position)
    elseif T<T_cr
        error("Crystallization Error! Temperature of LiBr in solution ("+T+"K) below min. value of "+T_cr+"K at position: "+position)
    end
end
end 