function NH3_pTDiagram(T, p, s)
% ----------------------------------------------------------------------- %
%{
Author  : Ludwig Irrgang
Date    : 01.07.2023
Copyright information:
Ludwig Irrgang
Lehrstuhl für Energiesysteme
TUM School of Engineering and Design
Technische Universität München
Boltzmannstr. 15 
85748 Garching b. München
ludwig.irrgang@tum.de
%}
% ------------------- Iso concentration lines --------------------------- %
% Define concentrations
w = [0;0.4;0.45;0.5;0.55;0.6;0.65;0.7;0.75;0.8;0.85;0.9;0.95;1];
% Define temperatures
T_vap_sol = 0:0.5:(T.sol_des_out-273.15);
if s.configuration == "Double_effect"
    T_vap_sol = 0:0.5:(T.sol_des_outI-273.15);
end
if s.configuration == "Double_lift"
    T_vap_sol = 0:0.5:(T.sol_des_outI-273.15);
end
T_vap_sol = T_vap_sol + 273.15;
% Initialize variables
p_vap_sol = zeros(length(w),length(T_vap_sol));
T_vap_plot = zeros(length(w),length(T_vap_sol));
% ------------------------- Cycle lines --------------------------------- %
switch s.configuration
    case 'Base'
        T_cond = CoolProp.PropsSI('T','P',p.cond,'Q',1,'Ammonia');
        T_evap = CoolProp.PropsSI('T','P',p.evap,'Q',1,'Ammonia');
        Cycle_plot_sol = [T.sol_abs_out,T.sol_des_in,T.sol_des_out,T.sol_abs_in,T.sol_abs_out;T_evap,T_cond,T_cond,T_evap,T_evap];
        Cycle_plot_ref = [T_cond,T_evap,T_cond];
        Connection_ref_high = [T_cond,T.sol_des_in;T_cond,T_cond];
        Connection_ref_low = [T_evap,T.sol_abs_out;T_evap,T_evap];
    case 'Double_effect'
        T_cond_int = CoolProp.PropsSI('T','P',p.cond_int,'Q',1,'Ammonia');
        T_cond = CoolProp.PropsSI('T','P',p.cond,'Q',1,'Ammonia');
        T_evap = CoolProp.PropsSI('T','P',p.evap,'Q',1,'Ammonia');
        Cycle_plot_sol = [T.sol_abs_out,T.sol_des_in,T.sol_des_out,T.sol_abs_in,T.sol_abs_out;T_evap,T_cond,T_cond,T_evap,T_evap];
        Cycle_plot_solI = [T.sol_abs_out,T.sol_des_inI,T.sol_des_outI,T.sol_abs_inI,T.sol_abs_out;T_evap,T_cond_int,T_cond_int,T_evap,T_evap];
        Cycle_plot_ref = [T_cond_int,T_evap,T_cond_int];
        Connection_ref_high = [T_cond_int,T.sol_des_inI;T_cond_int,T_cond_int];
        Connection_ref_mid = [T_cond,T.sol_des_in;T_cond,T_cond];
        Connection_ref_low = [T_evap,T.sol_abs_out;T_evap,T_evap];
    case 'Double_lift'
        T_cond = CoolProp.PropsSI('T','P',p.cond,'Q',1,'Ammonia');
        T_evap = CoolProp.PropsSI('T','P',p.evap,'Q',1,'Ammonia');
        T_mid = CoolProp.PropsSI('T','P',p.mid,'Q',1,'Ammonia');
        Cycle_plot_sol = [T.sol_abs_out,T.sol_des_in,T.sol_des_out,T.sol_abs_in,T.sol_abs_out;T_evap,T_mid,T_mid,T_evap,T_evap];
        Cycle_plot_solI = [T.sol_abs_outI,T.sol_des_inI,T.sol_des_outI,T.sol_abs_inI,T.sol_abs_outI;T_mid,T_cond,T_cond,T_mid,T_mid];
        Cycle_plot_ref = [T_cond,T_evap,T_cond];
        Connection_ref_high = [T_cond,T.sol_des_inI;T_cond,T_cond];
        Connection_ref_mid = [T.sol_abs_inI,T.sol_des_in;T_mid,T_mid];
        Connection_ref_low = [T_evap,T.sol_abs_out;T_evap,T_evap];
    case 'HP_comp_base'
        T_cond = CoolProp.PropsSI('T','P',p.cond,'Q',1,'Ammonia');
        T_des = CoolProp.PropsSI('T','P',p.des,'Q',1,'Ammonia');
        T_evap = CoolProp.PropsSI('T','P',p.evap,'Q',1,'Ammonia');
        Cycle_plot_sol = [T.sol_abs_out,T.sol_des_in,T.sol_des_out,T.sol_abs_in,T.sol_abs_out;T_evap,T_des,T_des,T_evap,T_evap];
        Cycle_plot_ref = [T_cond,T_evap,T_cond];
        Connection_ref_high = [T_cond,T.sol_des_in;T_cond,T_des];
        Connection_ref_low = [T_evap,T.sol_abs_out;T_evap,T_evap];
    case 'LP_comp_base'
        T_cond = CoolProp.PropsSI('T','P',p.cond,'Q',1,'Ammonia');
        T_abs = CoolProp.PropsSI('T','P',p.abs,'Q',1,'Ammonia');
        T_evap = CoolProp.PropsSI('T','P',p.evap,'Q',1,'Ammonia');
        Cycle_plot_sol = [T.sol_abs_out,T.sol_des_in,T.sol_des_out,T.sol_abs_in,T.sol_abs_out;T_abs,T_cond,T_cond,T_abs,T_abs];
        Cycle_plot_ref = [T_cond,T_evap,T_cond];
        Connection_ref_high = [T_cond,T.sol_des_in;T_cond,T_cond];
        Connection_ref_low = [T_evap,T.sol_abs_out;T_evap,T_abs];
    case 'Cascade_base_comp'
        if s.abs_active == 0
            return
        end
        T_cond = CoolProp.PropsSI('T','P',p.cond,'Q',1,'Ammonia');
        T_evap = CoolProp.PropsSI('T','P',p.evap,'Q',1,'Ammonia');
        Cycle_plot_sol = [T.sol_abs_out,T.sol_des_in,T.sol_des_out,T.sol_abs_in,T.sol_abs_out;T_evap,T_cond,T_cond,T_evap,T_evap];
        Cycle_plot_ref = [T_cond,T_evap,T_cond];
        Connection_ref_high = [T_cond,T.sol_des_in;T_cond,T_cond];
        Connection_ref_low = [T_evap,T.sol_abs_out;T_evap,T_evap];
    otherwise
        error("Configuration is not defined correctly")
end
% -------------------- Iso-concentration-lines -------------------------- %
for i=1:1:length(w)
    for j=1:1:length(T_vap_sol)
        try
            p_vap_sol(i,j) = refpropm('P','T',T_vap_sol(j),'Q',0,'AMMONIA','WATER',[1-w(i) w(i)])*1000;
            T_vap_plot(i,j) = CoolProp.PropsSI('T','P',p_vap_sol(i,j),'Q',1,'Ammonia');
        catch
            p_vap_sol(i,j) = nan;
            T_vap_plot(i,j) = nan;
        end
        if T_vap_plot(i,j) < 273.15
            p_vap_sol(i,j) = nan;
            T_vap_plot(i,j) = nan;
        end
    end
end
% save('10_Plots/Rev_Duhring_LiBr.mat');
% ------------------------------ Plot ----------------------------------- %
figure
plot(T_vap_sol(:)-273.15,T_vap_plot(find(w==0),:)-273.15,'b','DisplayName',"Refrigerant",'LineWidth',2.0)
hold on
plot(T_vap_sol(:)-273.15,T_vap_plot(find(w==0.4),:)-273.15,'DisplayName',"w_{H_2O}=0,4",'LineWidth',1.0)
hold on
plot(T_vap_sol(:)-273.15,T_vap_plot(find(w==0.45),:)-273.15,'k','HandleVisibility','off','LineWidth',0.01)
hold on
plot(T_vap_sol(:)-273.15,T_vap_plot(find(w==0.5),:)-273.15,'k','HandleVisibility','off','LineWidth',0.01)
hold on
plot(T_vap_sol(:)-273.15,T_vap_plot(find(w==0.55),:)-273.15,'k','HandleVisibility','off','LineWidth',0.01)
hold on
plot(T_vap_sol(:)-273.15,T_vap_plot(find(w==0.6),:)-273.15,'DisplayName',"w_{H_2O}=0,6",'LineWidth',1.0)
hold on
plot(T_vap_sol(:)-273.15,T_vap_plot(find(w==0.65),:)-273.15,'k','HandleVisibility','off','LineWidth',0.01)
hold on
plot(T_vap_sol(:)-273.15,T_vap_plot(find(w==0.7),:)-273.15,'k','HandleVisibility','off','LineWidth',0.01)
hold on
plot(T_vap_sol(:)-273.15,T_vap_plot(find(w==0.75),:)-273.15,'k','HandleVisibility','off','LineWidth',0.01)
hold on
plot(T_vap_sol(:)-273.15,T_vap_plot(find(w==0.8),:)-273.15,'DisplayName',"w_{H_2O}=0,8",'LineWidth',1.0)
hold on
plot(T_vap_sol(:)-273.15,T_vap_plot(find(w==0.85),:)-273.15,'k','HandleVisibility','off','LineWidth',0.01)
hold on
plot(T_vap_sol(:)-273.15,T_vap_plot(find(w==0.9),:)-273.15,'k','HandleVisibility','off','LineWidth',0.01)
hold on
plot(T_vap_sol(:)-273.15,T_vap_plot(find(w==0.95),:)-273.15,'k','HandleVisibility','off','LineWidth',0.01)
hold on
plot(T_vap_sol(:)-273.15,T_vap_plot(find(w==1),:)-273.15,'DisplayName',"w_{H_2O}=1",'LineWidth',1.0)
hold on
% ----------------------------- cycle ----------------------------------- %
switch s.configuration
    case 'Base'
        plot(Cycle_plot_sol(1,:)-273.15,Cycle_plot_sol(2,:)-273.15,'m','DisplayName',"Solution Cycle",'LineWidth',3.0)
        hold on
        plot(Cycle_plot_ref-273.15,Cycle_plot_ref-273.15,'g','DisplayName',"Refrigerant Cycle",'LineWidth',3.0)
        hold on
        plot(Connection_ref_high(1,:)-273.15,Connection_ref_high(2,:)-273.15,'k--','HandleVisibility','off','LineWidth',1.0)
        hold on
        plot(Connection_ref_low(1,:)-273.15,Connection_ref_low(2,:)-273.15,'k--','HandleVisibility','off','LineWidth',1.0)
    case 'Double_effect'
        plot(Cycle_plot_sol(1,:)-273.15,Cycle_plot_sol(2,:)-273.15,'m','DisplayName',"Solution Cycle",'LineWidth',3.0)
        hold on
        plot(Cycle_plot_solI(1,:)-273.15,Cycle_plot_solI(2,:)-273.15,'m','HandleVisibility','off','LineWidth',3.0)
        hold on
        plot(Cycle_plot_ref-273.15,Cycle_plot_ref-273.15,'g','DisplayName',"Refrigerant Cycle",'LineWidth',3.0)
        hold on
        plot(Connection_ref_high(1,:)-273.15,Connection_ref_high(2,:)-273.15,'k--','HandleVisibility','off','LineWidth',1.0)
        hold on
        plot(Connection_ref_mid(1,:)-273.15,Connection_ref_mid(2,:)-273.15,'k--','HandleVisibility','off','LineWidth',1.0)
        hold on
        plot(Connection_ref_low(1,:)-273.15,Connection_ref_low(2,:)-273.15,'k--','HandleVisibility','off','LineWidth',1.0)
    case 'Double_lift'
        plot(Cycle_plot_sol(1,:)-273.15,Cycle_plot_sol(2,:)-273.15,'m','DisplayName',"Solution Cycle",'LineWidth',3.0)
        hold on
        plot(Cycle_plot_solI(1,:)-273.15,Cycle_plot_solI(2,:)-273.15,'m','HandleVisibility','off','LineWidth',3.0)
        hold on
        plot(Cycle_plot_ref-273.15,Cycle_plot_ref-273.15,'g','DisplayName',"Refrigerant Cycle",'LineWidth',3.0)
        hold on
        plot(Connection_ref_high(1,:)-273.15,Connection_ref_high(2,:)-273.15,'k--','HandleVisibility','off','LineWidth',1.0)
        hold on
        plot(Connection_ref_mid(1,:)-273.15,Connection_ref_mid(2,:)-273.15,'k--','HandleVisibility','off','LineWidth',1.0)
        hold on
        plot(Connection_ref_low(1,:)-273.15,Connection_ref_low(2,:)-273.15,'k--','HandleVisibility','off','LineWidth',1.0)
    case 'HP_comp_base'
        plot(Cycle_plot_sol(1,:)-273.15,Cycle_plot_sol(2,:)-273.15,'m','DisplayName',"Solution Cycle",'LineWidth',3.0)
        hold on
        plot(Cycle_plot_ref-273.15,Cycle_plot_ref-273.15,'g','DisplayName',"Refrigerant Cycle",'LineWidth',3.0)
        hold on
        plot(Connection_ref_high(1,:)-273.15,Connection_ref_high(2,:)-273.15,'k--','HandleVisibility','off','LineWidth',1.0)
        hold on
        plot(Connection_ref_low(1,:)-273.15,Connection_ref_low(2,:)-273.15,'k--','HandleVisibility','off','LineWidth',1.0)
    case 'LP_comp_base'
        plot(Cycle_plot_sol(1,:)-273.15,Cycle_plot_sol(2,:)-273.15,'m','DisplayName',"Solution Cycle",'LineWidth',3.0)
        hold on
        plot(Cycle_plot_ref-273.15,Cycle_plot_ref-273.15,'g','DisplayName',"Refrigerant Cycle",'LineWidth',3.0)
        hold on
        plot(Connection_ref_high(1,:)-273.15,Connection_ref_high(2,:)-273.15,'k--','HandleVisibility','off','LineWidth',1.0)
        hold on
        plot(Connection_ref_low(1,:)-273.15,Connection_ref_low(2,:)-273.15,'k--','HandleVisibility','off','LineWidth',1.0)
    case 'Cascade_base_comp'
        plot(Cycle_plot_sol(1,:)-273.15,Cycle_plot_sol(2,:)-273.15,'m','DisplayName',"Solution Cycle",'LineWidth',3.0)
        hold on
        plot(Cycle_plot_ref-273.15,Cycle_plot_ref-273.15,'g','DisplayName',"Refrigerant Cycle",'LineWidth',3.0)
        hold on
        plot(Connection_ref_high(1,:)-273.15,Connection_ref_high(2,:)-273.15,'k--','HandleVisibility','off','LineWidth',1.0)
        hold on
        plot(Connection_ref_low(1,:)-273.15,Connection_ref_low(2,:)-273.15,'k--','HandleVisibility','off','LineWidth',1.0)
end
legend('Location','northwest');
xlabel('T_{sol,boil}, [°C]');
ylabel('T_{ref,boil}, [°C]');
hold off
title('Cycle Temperatures in Concentration Daigram');
grid minor
% ----------------------------------------------------------------------- %
end