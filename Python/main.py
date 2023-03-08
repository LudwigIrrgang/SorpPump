# ----------------------------------------------------------------------- #
# Import the relevant modules and libraries
# ----------------------------------------------------------------------- #
import time
start_time = time.time()
import sys
import numpy 
import pandas
import math
from scipy.optimize import fsolve, newton
import CoolProp.CoolProp 
# ----------------------------------------------------------------------- #
# Import the Fluids and Internal_Cycle_model modules
# ----------------------------------------------------------------------- #
from Fluids import LiBrSol, H2O 
from Internal_Cycle_models import base_model_H2OLiBr, base_model_NH3H2O, doubleEffect_model_H2OLiBr, doubleLift_model_H2OLiBr, doubleLift_model_NH3H2O, doubleEffect_model_NH3H2O
print("--- %s seconds ---" % (time.time() - start_time))

# ----------------------------------------------------------------------- #
# Define functions for printing and storing the results
# ----------------------------------------------------------------------- #

def print_results_to_console(T, p, h, m, w, eta, Q, PP, s):
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
    print(f"----T:----\n{df_T}\n ----p:----\n{df_p}\n ----h:----\n{df_h} \n ----m:----\n{df_m} \n ----w:----\n{df_w}")


def save_results_to_csv(T, p, h, m, w, eta, Q, PP, s):
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
    string_csv = input("name or location for storing the result .csv file without \".csv\": ")
    result_folder = "Python\Results"
    df_result.to_csv(result_folder + "\\" + string_csv +".csv", sep=";")
    print('The results are stored in: ' + result_folder)

# ----------------------------------------------------------------------- #
# Define the class for handeling the variables, and initialise them
# ----------------------------------------------------------------------- #
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

# ----------------------------------------------------------------------- #
# Make input decission: manual in python-file or step py step through input()
# ----------------------------------------------------------------------- #

external_input = True   # if set to False all information have to be changed in the python file

# ----------------------------------------------------------------------- #
# ----------------------------------------------------------------------- #
# Make the inputs for the model
# ----------------------------------------------------------------------- #
# ----------------------------------------------------------------------- #
if external_input ==False:

    working_fluid =    "NH3_H2O"    # Choose the working fluid ("NH3_H2O" or "LiBr_H2O")

    cycle_model =      "DLK"       # Choose the cycle_model ("base", "DEK" or "DLK")

    save_as_csv =      "True"       # Save the results as .csv? ("True" or "False")

    # Temperature informations needed:
    T.evap= 278.15              
    T.sol_abs_out= 306.15      
    T.sol_des_out = 349.77  
    T.cond = 309.15          
    T.ext_cond_in = 303.15 
    T.cond_int =   352.77
    

    # Defining which heat flow will be provided. Q_evap or Q_des?

    s.requirement = "Q_evap"

    # Value of heat flow
    Q.dec = 10_000             # Q_des

    # Efficiency of the pump
    eta.pump = 1    

    #Heat exchanger parameters
    HX.T_PP_SHEX = 3
    HX.T_PP_RHEX = 3
    HX.T_PP_cond = 3
    HX.T_PP_cond = 3
    HX.T_PP_SHEXI =3
    HX.T_PP_cond_int = 3

    ## Additional Parameters needed based on the differnt cycle types an working fluids
    # For both LiBr-H2O and NH3-H2O

    # Based on cycle, see discription behind value

    T.cond_int =   352.77       # DEK
    T.sol_des_outI= 410.15      # DEK & DLK
    T.sol_abs_outI = 306.15     # DLK 

    HX.T_PP_SHEXI = 3           # DEK & DLK
    HX.T_PP_cond_int = 3        # DEK

    ## only for NH3-H2O
    # all cycle types
    HX.dT_ref_des = 5

    # Double Effect + Double Lift
    HX.dT_ref_desI = 5
      
elif external_input == True:
    # ----------------------------------------------------------------------- #
    # requesting all needed inputs 
    # ----------------------------------------------------------------------- #
    working_fluid = input("Choose the working fluid (\"NH3_H2O\" or \"LiBr_H2O\") :")

    cycle_model = input("Choose the cycle_model (\"base\", \"DEK\" or \"DLK\") :")

    save_as_csv = input("Save the results as .csv? (True or False):")

    # ----------------------------------------------------------------------- #
    # requesting inputs for input data
    # ----------------------------------------------------------------------- #
    # input data which is needed for all cycle types

    T.evap = float(input("T.evap = ? [K]"))
    T.sol_abs_out = float(input("T.sol_abs_out = ? [K]"))
    T.sol_des_out = float(input("T.sol_des_out = ? [K]"))
    T.cond = float(input("T.cond = ? [K]"))
    T.ext_cond_in = float(input("T.ext_cond_in [K]"))

    s.requirement = input("Which heat flow will be provided? (Q_evap or Q_des)")

    Q.dec = float(input("Q_evap or Q_des based on s.requirement = ? [W]"))

    eta.pump = float(input("Efficiecy of the pump = ? "))

    HX.T_PP_SHEX = float(input("Heat exchanger parameter HX.T_PP_SHEX = ? "))
    HX.T_PP_RHEX = float(input("Heat exchanger parameter HX.T_PP_RHEX = ? "))
    HX.T_PP_cond = float(input("Heat exchanger parameter HX.T_PP_cond = ? "))

    if cycle_model == "DEK":
        T.sol_des_outI = float(input("T.sol_des_outI = ? [K]"))

        HX.T_PP_SHEXI = float(input("Heat exchanger parameter HX.T_PP_SHEXI = ? "))
        HX.T_PP_cond_int = float(input("Heat exchanger parameter HX.T_PP_cond_int = ? "))

        if working_fluid =="NH3_H2O":
            HX.dT_ref_des = float(input("Heat exchanger parameter HX.dT_ref_des= ? "))
            HX.dT_ref_desI = float(input("Heat exchanger parameter HX.dT_ref_desI= ? "))
    
    elif cycle_model == "DLK":
        T.sol_abs_outI = float(input("T.sol_abs_outI = ? [K]"))
        T.sol_des_outI = float(input("T.sol_des_outI = ? [K]"))

        HX.T_PP_SHEXI = float(input("Heat exchanger parameter HX.T_PP_SHEXI = ? "))

        if working_fluid =="NH3_H2O":
            HX.dT_ref_des = float(input("Heat exchanger parameter HX.dT_ref_des= ? "))
            HX.dT_ref_desI = float(input("Heat exchanger parameter HX.dT_ref_desI= ? "))

    elif cycle_model == "base":
        if working_fluid =="NH3_H2O":
            HX.dT_ref_des = float(input("Heat exchanger parameter HX.dT_ref_des= ? "))
    else:
        sys.exit("Cycle model ist not defined!")  


if cycle_model == "base":
    
    if working_fluid == "LiBr_H2O":
        T, p, h, m, w, eta, Q, PP, s= base_model_H2OLiBr.base_model_H2OLiBr(T, p, h, m, eta, Q, HX, s)
        print("------------------------------base_model_H2O_LiBr------------------------------")

        print_results_to_console(T, p, h, m, w, eta, Q, PP, s)

        if save_as_csv == "True":
            save_results_to_csv(T, p, h, m, w, eta, Q, PP, s)
    
    elif working_fluid == "NH3_H2O":
        T, p, h, m, w, eta, Q, PP, s= base_model_NH3H2O.base_model_NH3H2O(T, p, h, m, eta, Q, HX, s)
        print("------------------------------base_model_NH3_H2O------------------------------")

        print_results_to_console(T, p, h, m, w, eta, Q, PP, s)

        if save_as_csv == "True":
            save_results_to_csv(T, p, h, m, w, eta, Q, PP, s)
    
    else:
        sys.exit("Working fluid is not defined!")

elif cycle_model =="DEK":
    if working_fluid == "LiBr_H2O":
        T, p, h, m, w, eta, Q, PP, s= doubleEffect_model_H2OLiBr.doubleEffect_model_H2OLiBr(T, p, h, m, eta, Q, HX, s)
        print("------------------------------double_effect_H2O_LiBr------------------------------")

        print_results_to_console(T, p, h, m, w, eta, Q, PP, s)

        if save_as_csv == "True":
            save_results_to_csv(T, p, h, m, w, eta, Q, PP, s)
    
    elif working_fluid == "NH3_H2O":
        T, p, h, m, w, eta, Q, PP, s= doubleEffect_model_NH3H2O.doubleEffect_model_NH3H2O(T, p, h, m, eta, Q, HX, s)
        print("------------------------------double_effect_NH3_H2O------------------------------")

        print_results_to_console(T, p, h, m, w, eta, Q, PP, s)

        if save_as_csv == "True":
            save_results_to_csv(T, p, h, m, w, eta, Q, PP, s)
    
    else:
        sys.exit("Working fluid is not defined!")

elif cycle_model =="DLK":
    if working_fluid == "LiBr_H2O":
        T, p, h, m, w, eta, Q, PP, s= doubleLift_model_H2OLiBr.doubleLift_model_H2OLiBr(T, p, h, m, eta, Q, HX, s)
        print("------------------------------double_lift_H2O_LiBr------------------------------")

        print_results_to_console(T, p, h, m, w, eta, Q, PP, s)

        if save_as_csv == "True":
            save_results_to_csv(T, p, h, m, w, eta, Q, PP, s)
    
    elif working_fluid == "NH3_H2O":
        T, p, h, m, w, eta, Q, PP, s= doubleLift_model_NH3H2O.doubleLift_model_NH3H2O(T, p, h, m, eta, Q, HX, s)
        print("------------------------------double_lift_NH3_H2O------------------------------")

        print_results_to_console(T, p, h, m, w, eta, Q, PP, s)

        if save_as_csv == "True":
            save_results_to_csv(T, p, h, m, w, eta, Q, PP, s)
    
    else:
        sys.exit("Working fluid is not defined!")


        




