import time
start_time = time.time()
import numpy 
import pandas
import math
from scipy.optimize import fsolve, newton
import CoolProp.CoolProp 
from Fluids import LiBrSol, H2O 
from Internal_Cycle_models import base_model_H2OLiBr, base_model_NH3H2O, doubleEffect_model_H2OLiBr, doubleLift_model_H2OLiBr
print("--- %s seconds ---" % (time.time() - start_time))

class var():
    i =0

T   =   var()
p   =   var()
h   =   var()
m   =   var()
eta =   var()
Q   =   var()
HX  =   var()
s   =   var()

# initialising the requiered input data
T.evap= 278.15              
T.sol_abs_out= 306.15      
T.sol_des_out = 350.15    
T.cond = 309.15          
T.ext_cond_in = 303.15     

Q.dec = 10_000             # Q_des

eta.pump = 1    

HX.T_PP_SHEX = 3
HX.T_PP_RHEX = 3
HX.T_PP_cond = 3

s.requirement = "Q_evap"





T, p, h, m, w, eta, Q, PP, s= base_model_H2OLiBr.base_model_H2OLiBr(T, p, h, m, eta, Q, HX, s)
print("------------------------------base_model_H2OLiBr------------------------------")
print("--- %s seconds ---" % (time.time() - start_time))
df_T = pandas.DataFrame.from_dict(T.__dict__,orient="index")
df_T["propertie"]="T"
df_p = pandas.DataFrame.from_dict(p.__dict__,orient="index")
df_p["propertie"]="p"
df_h = pandas.DataFrame.from_dict(h.__dict__,orient="index")
df_h["propertie"] = "h"
df_m = pandas.DataFrame.from_dict(m.__dict__,orient="index")
df_m["propertie"] = "m"
df_w = pandas.DataFrame.from_dict(w.__dict__,orient="index")
df_w["propertie"] ="w"
df_eta = pandas.DataFrame.from_dict(eta.__dict__,orient="index")
df_eta["propertie"] = "eta"
df_Q = pandas.DataFrame.from_dict(Q.__dict__,orient="index")
df_Q["propertie"] = "Q"
df_PP = pandas.DataFrame.from_dict(PP.__dict__,orient="index")
df_PP["propertie"] = "PP"
df_s = pandas.DataFrame.from_dict(s.__dict__,orient="index")
df_s["propertie"] = "s"

df_result = pandas.concat([df_T,df_p,df_h,df_m,df_w, df_eta, df_Q, df_PP, df_s ])
df_result.to_csv("result_LiBr_H2O_base_model_test.csv")
print(f"----T:----\n{df_T}\n ----p:----\n{df_p}\n ----h:----\n{df_h} \n ----m:----\n{df_m} \n ----w:----\n{df_w}")
print(f"------\n------\n---LiBr_doubleEffect---\n")
T   =   var()
p   =   var()
h   =   var()
m   =   var()
eta =   var()
Q   =   var()
HX  =   var()
s   =   var()

# initialising the requiered input data

T.evap= 278.15              
T.sol_abs_out= 306.15      
T.sol_des_out = 350.15    
T.cond = 309.15          
T.ext_cond_in = 303.15 
T.cond_int =   350.4934
T.sol_des_outI=410.15 

Q.dec = 10_000             # Q_des

eta.pump = 1    

HX.T_PP_SHEX = 3
HX.T_PP_RHEX = 3
HX.T_PP_cond = 3
HX.T_PP_cond = 3
HX.T_PP_SHEXI =3
HX.T_PP_cond_int = 3

s.requirement = "Q_evap"





T, p, h, m, w, eta, Q, PP, s= doubleEffect_model_H2OLiBr.doubleEffect_model_H2OLiBr(T, p, h, m, eta, Q, HX, s)
print("------------------------------double Effect_H2OLiBr------------------------------")
print("--- %s seconds ---" % (time.time() - start_time))
df_T = pandas.DataFrame.from_dict(T.__dict__,orient="index")
df_T["propertie"]="T"
df_p = pandas.DataFrame.from_dict(p.__dict__,orient="index")
df_p["propertie"]="p"
df_h = pandas.DataFrame.from_dict(h.__dict__,orient="index")
df_h["propertie"] = "h"
df_m = pandas.DataFrame.from_dict(m.__dict__,orient="index")
df_m["propertie"] = "m"
df_w = pandas.DataFrame.from_dict(w.__dict__,orient="index")
df_w["propertie"] ="w"
df_eta = pandas.DataFrame.from_dict(eta.__dict__,orient="index")
df_eta["propertie"] = "eta"
df_Q = pandas.DataFrame.from_dict(Q.__dict__,orient="index")
df_Q["propertie"] = "Q"
df_PP = pandas.DataFrame.from_dict(PP.__dict__,orient="index")
df_PP["propertie"] = "PP"
df_s = pandas.DataFrame.from_dict(s.__dict__,orient="index")
df_s["propertie"] = "s"

df_result = pandas.concat([df_T,df_p,df_h,df_m,df_w, df_eta, df_Q, df_PP, df_s ])
df_result.to_csv("result_LiBr_H2O_doubleEffect_model_test.csv")
print(f"----T:----\n{df_T}\n ----p:----\n{df_p}\n ----h:----\n{df_h} \n ----m:----\n{df_m} \n ----w:----\n{df_w}")

print(f"------\n------\n---LiBr_doubleLift---\n")
T   =   var()
p   =   var()
h   =   var()
m   =   var()
eta =   var()
Q   =   var()
HX  =   var()
s   =   var()

# initialising the requiered input data

T.evap= 278.15              
T.sol_abs_out= 306.15      
T.sol_des_out = 350.15    
T.cond = 309.15          
T.ext_cond_in = 303.15 

T.sol_abs_outI = 306.15
T.sol_des_outI= 330.15


Q.dec = 10_000             # Q_des

eta.pump = 1    

HX.T_PP_SHEX = 3
HX.T_PP_RHEX = 3
HX.T_PP_cond = 3
HX.T_PP_SHEXI =3


s.requirement = "Q_evap"





T, p, h, m, w, eta, Q, PP, s= doubleLift_model_H2OLiBr.doubleLift_model_H2OLiBr(T, p, h, m, eta, Q, HX, s)
print("------------------------------double Lift_H2OLiBr------------------------------")
print("--- %s seconds ---" % (time.time() - start_time))
df_T = pandas.DataFrame.from_dict(T.__dict__,orient="index")
df_T["propertie"]="T"
df_p = pandas.DataFrame.from_dict(p.__dict__,orient="index")
df_p["propertie"]="p"
df_h = pandas.DataFrame.from_dict(h.__dict__,orient="index")
df_h["propertie"] = "h"
df_m = pandas.DataFrame.from_dict(m.__dict__,orient="index")
df_m["propertie"] = "m"
df_w = pandas.DataFrame.from_dict(w.__dict__,orient="index")
df_w["propertie"] ="w"
df_eta = pandas.DataFrame.from_dict(eta.__dict__,orient="index")
df_eta["propertie"] = "eta"
df_Q = pandas.DataFrame.from_dict(Q.__dict__,orient="index")
df_Q["propertie"] = "Q"
df_PP = pandas.DataFrame.from_dict(PP.__dict__,orient="index")
df_PP["propertie"] = "PP"
df_s = pandas.DataFrame.from_dict(s.__dict__,orient="index")
df_s["propertie"] = "s"

df_result = pandas.concat([df_T,df_p,df_h,df_m,df_w, df_eta, df_Q, df_PP, df_s ])
df_result.to_csv("result_LiBr_H2O_doubleLift_model_test.csv")
print(f"----T:----\n{df_T}\n ----p:----\n{df_p}\n ----h:----\n{df_h} \n ----m:----\n{df_m} \n ----w:----\n{df_w}")


print(f"------\n------\n---NH3_H2O_base_model---\n")
T   =   var()
p   =   var()
h   =   var()
m   =   var()
eta =   var()
Q   =   var()
HX  =   var()
s   =   var()

# initialising the requiered input data

T.evap= 278.15              
T.sol_abs_out= 306.15      
T.sol_des_out = 350.15    
T.cond = 309.15          
T.ext_cond_in = 303.15 

Q.dec = 10_000             # Q_des

eta.pump = 1    

HX.T_PP_SHEX = 3
HX.T_PP_RHEX = 3
HX.T_PP_cond = 3

s.requirement = "Q_evap"


T, p, h, m, w, eta, Q, PP, s= base_model_NH3H2O.base_model_NH3H2O(T, p, h, m, eta, Q, HX, s)
print("------------------------------base_NH3H2O------------------------------")
print("--- %s seconds ---" % (time.time() - start_time))
df_T = pandas.DataFrame.from_dict(T.__dict__,orient="index")
df_T["propertie"]="T"
df_p = pandas.DataFrame.from_dict(p.__dict__,orient="index")
df_p["propertie"]="p"
df_h = pandas.DataFrame.from_dict(h.__dict__,orient="index")
df_h["propertie"] = "h"
df_m = pandas.DataFrame.from_dict(m.__dict__,orient="index")
df_m["propertie"] = "m"
df_w = pandas.DataFrame.from_dict(w.__dict__,orient="index")
df_w["propertie"] ="w"
df_eta = pandas.DataFrame.from_dict(eta.__dict__,orient="index")
df_eta["propertie"] = "eta"
df_Q = pandas.DataFrame.from_dict(Q.__dict__,orient="index")
df_Q["propertie"] = "Q"
df_PP = pandas.DataFrame.from_dict(PP.__dict__,orient="index")
df_PP["propertie"] = "PP"
df_s = pandas.DataFrame.from_dict(s.__dict__,orient="index")
df_s["propertie"] = "s"

df_result = pandas.concat([df_T,df_p,df_h,df_m,df_w, df_eta, df_Q, df_PP, df_s ])
df_result.to_csv("result_NH3_H2O_base_model_test.csv")
print(f"----T:----\n{df_T}\n ----p:----\n{df_p}\n ----h:----\n{df_h} \n ----m:----\n{df_m} \n ----w:----\n{df_w}")
