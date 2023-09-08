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
# ----------------------------------------------------------------------- #
# Define functions for printing and storing the results
# ----------------------------------------------------------------------- #


def print_results_to_console(T, p, h, m, w, eta, Q, PP, s):
    df_T = pandas.DataFrame.from_dict(T.__dict__, orient="index")
    df_T["propertie"]="T"
    df_p = pandas.DataFrame.from_dict(p.__dict__, orient="index")
    df_p["propertie"]="p"
    df_h = pandas.DataFrame.from_dict(h.__dict__, orient="index")
    df_h["propertie"] = "h"
    df_m = pandas.DataFrame.from_dict(m.__dict__, orient="index")
    df_m["propertie"] = "m"
    df_w = pandas.DataFrame.from_dict(w.__dict__, orient="index")
    df_w["propertie"] ="w"
    df_eta = pandas.DataFrame.from_dict(eta.__dict__, orient="index")
    df_eta["propertie"] = "eta"
    df_Q = pandas.DataFrame.from_dict(Q.__dict__, orient="index")
    df_Q["propertie"] = "Q"
    df_PP = pandas.DataFrame.from_dict(PP.__dict__, orient="index")
    df_PP["propertie"] = "PP"
    df_s = pandas.DataFrame.from_dict(s.__dict__, orient="index")
    df_s["propertie"] = "s"

    df_result = pandas.concat([df_T, df_p, df_h, df_m, df_w, df_eta, df_Q, df_PP, df_s ])
    print(f"----T:----\n{df_T}\n ----p:----\n{df_p}\n ----h:----\n{df_h} \n ----m:----\n{df_m} \n ----w:----\n{df_w} \n ----Q:----\n{df_Q} \n ----Q:----\n{df_PP}")


def save_results_to_csv(T, p, h, m, w, eta, Q, PP, s):
    df_T = pandas.DataFrame.from_dict(T.__dict__, orient="index")
    df_T["propertie"]="T"
    df_p = pandas.DataFrame.from_dict(p.__dict__, orient="index")
    df_p["propertie"]="p"
    df_h = pandas.DataFrame.from_dict(h.__dict__, orient="index")
    df_h["propertie"] = "h"
    df_m = pandas.DataFrame.from_dict(m.__dict__, orient="index")
    df_m["propertie"] = "m"
    df_w = pandas.DataFrame.from_dict(w.__dict__, orient="index")
    df_w["propertie"] ="w"
    df_eta = pandas.DataFrame.from_dict(eta.__dict__, orient="index")
    df_eta["propertie"] = "eta"
    df_Q = pandas.DataFrame.from_dict(Q.__dict__, orient="index")
    df_Q["propertie"] = "Q"
    df_PP = pandas.DataFrame.from_dict(PP.__dict__, orient="index")
    df_PP["propertie"] = "PP"
    df_s = pandas.DataFrame.from_dict(s.__dict__, orient="index")
    df_s["propertie"] = "s"
    df_result = pandas.concat([df_T, df_p, df_h, df_m, df_w, df_eta, df_Q, df_PP, df_s])
    string_csv = input("name or location for storing the result .csv file without \".csv\": ")
    result_folder = "Results"
    df_result.to_csv(result_folder + "\\" + string_csv + ".csv", sep=";")
    print('The results are stored in: ' + result_folder)
# ----------------------------------------------------------------------- #
# Define the class for handling the variables, and initialise them
# ----------------------------------------------------------------------- #
class var():
    i = 0
T = var()
p = var()
h = var()
m = var()
eta = var()
Q = var()
HX = var()
s = var()
# ----------------------------------------------------------------------- #
# Make input decision: manual in python-file or step py step through input()
# ----------------------------------------------------------------------- #
external_input = False   # False - Input in file; True - Input through Console
# ----------------------------------------------------------------------- #
# ----------------------------------------------------------------------- #
# Make the inputs for the model
# ----------------------------------------------------------------------- #
# ----------------------------------------------------------------------- #
if not external_input:

    working_fluid = "LiBr_H2O"  # Choose the working fluid ("NH3_H2O" or "LiBr_H2O")

    cycle_model = "DE"        # Choose the cycle_model ("base", "DE" or "DL")

    save_as_csv = "True"        # Save the results as .csv? ("True" or "False")

    # External temperature boundary conditions
    # Heat source temperature
    T.ext_des_in = 150 + 273.15
    # Cold output temperature
    T.ext_evap_out = 6 + 273.15
    # Heat sink temperature
    T.ext_abs_in = 30 + 273.15
    T.ext_cond_in = 30 + 273.15

    # Heat exchanger PP temperature differences
    HX.T_PP_evap = 5
    HX.T_PP_abs = 3
    HX.T_PP_des = 3
    HX.T_PP_cond = 3
    HX.SC_cond = 0          # Defines subcooling at condenser
    HX.T_PP_SHEX = 5
    HX.T_PP_RHEX = 100
    HX.T_PP_SHEXI = 3       # Only necessary for DE and DL
    HX.T_PP_cond_int = 3    # Only necessary for DE
    HX.dT_ref_des = 5       # Defines superheating of refrigerant after desorber
    HX.dT_ref_desI = 5      # Defines superheating of refrigerant after desorber I

    # Defining which heat flow will be provided. Q_evap or Q_des?
    s.requirement = "Q_evap"

    # Value of heat flow
    Q.dec = 10000

    # Efficiency of the pump
    eta.pump = 1


elif external_input:
    # ----------------------------------------------------------------------- #
    # requesting all needed inputs 
    # ----------------------------------------------------------------------- #
    working_fluid = input("Choose the working fluid (\"NH3_H2O\" or \"LiBr_H2O\") :")

    cycle_model = input("Choose the cycle_model (\"base\", \"DE\" or \"DL\") :")

    save_as_csv = input("Save the results as .csv? (True or False):")

    # ----------------------------------------------------------------------- #
    # requesting inputs for input data
    # ----------------------------------------------------------------------- #
    # Input data which is needed for all cycle types

    T.ext_evap_out = float(input("Cold output = ? [K]"))
    T.ext_abs_in = float(input("Heat sink temperature = ? [K]"))
    T.ext_cond_in = T.ext_abs_in
    T.ext_des_in = float(input("Heat source temperature = ? [K]"))

    s.requirement = input("Which heat flow will be provided? (Q_evap or Q_des)")

    Q.dec = float(input("Q_evap or Q_des based on s.requirement = ? [W]"))

    eta.pump = float(input("Efficiency of the pump = ? "))

    HX.T_PP_SHEX = float(input("Heat exchanger parameter HX.T_PP_SHEX = ? "))
    HX.T_PP_RHEX = float(input("Heat exchanger parameter HX.T_PP_RHEX = ? "))
    HX.T_PP_cond = float(input("Heat exchanger parameter HX.T_PP_cond = ? "))

    HX.T_PP_evap = float(input("Evaporator PP temperature difference = ? [K]"))
    HX.T_PP_abs = float(input("Absorber PP temperature difference = ? [K]"))
    HX.T_PP_des = float(input("Desorber PP temperature difference = ? [K]"))
    HX.T_PP_cond = float(input("Condenser PP temperature difference = ? [K]"))
    HX.SC_cond = float(input("Sub-cooling at condenser = ? [K]"))
    HX.T_PP_SHEX = float(input("SHEX PP temperature difference = ? [K]"))
    HX.T_PP_RHEX = float(input("RHEX PP temperature difference = ? [K]"))
    HX.dT_ref_des = float(input("Refrigerant - solution temperature difference leaving desorber = ? [K]"))

    if cycle_model == "DE":
        HX.T_PP_SHEXI = float(input("SHEXI PP temperature difference = ? [K]"))
        HX.T_PP_cond_int = float(input("Internal condenser PP temperature difference = ? [K]"))
        HX.dT_ref_desI = float(input("Refrigerant - solution temperature difference leaving desorber I = ? [K]"))

    elif cycle_model == "DL":
        HX.T_PP_SHEXI = float(input("SHEXI PP temperature difference = ? [K]"))
        HX.dT_ref_desI = float(input("Refrigerant - solution temperature difference leaving desorber I = ? [K]"))

    else:
        sys.exit("Cycle model ist not defined!")


# Calculate Internal temperatures for cycle model input
T.sol_des_out = T.ext_des_in - HX.T_PP_SHEX
T.evap = T.ext_evap_out - HX.T_PP_evap
T.cond = T.ext_cond_in - HX.T_PP_cond - HX.SC_cond
T.sol_abs_out = T.ext_abs_in + HX.T_PP_abs

if cycle_model == "DE":
    # Internal condenser
    T.cond_int = T.cond + HX.dT_ref_des + HX.T_PP_cond_int
    # Upper pressure desorber
    T.sol_des_outI = T.sol_des_out
    # Middle pressure desorber
    T.sol_des_out = T.cond_int - HX.T_PP_cond_int
elif cycle_model == "DL":
    # Upper pressure desorber
    T.sol_des_outI = T.sol_des_out
    # Middle pressure absorber
    T.sol_abs_outI = T.sol_abs_out
else:
    sys.exit("Cycle model ist not defined!")


if cycle_model == "base":
    
    if working_fluid == "LiBr_H2O":
        T, p, h, m, w, eta, Q, PP, s = base_model_H2OLiBr.base_model_H2OLiBr(T, p, h, m, eta, Q, HX, s)
        print("------------------------------base_model_H2O_LiBr------------------------------")

        print_results_to_console(T, p, h, m, w, eta, Q, PP, s)

        if save_as_csv == "True":
            save_results_to_csv(T, p, h, m, w, eta, Q, PP, s)
    
    elif working_fluid == "NH3_H2O":
        T, p, h, m, w, eta, Q, PP, s = base_model_NH3H2O.base_model_NH3H2O(T, p, h, m, eta, Q, HX, s)
        print("------------------------------base_model_NH3_H2O------------------------------")

        print_results_to_console(T, p, h, m, w, eta, Q, PP, s)

        if save_as_csv == "True":
            save_results_to_csv(T, p, h, m, w, eta, Q, PP, s)
    
    else:
        sys.exit("Working fluid is not defined!")

elif cycle_model == "DE":
    if working_fluid == "LiBr_H2O":
        T, p, h, m, w, eta, Q, PP, s = doubleEffect_model_H2OLiBr.doubleEffect_model_H2OLiBr(T, p, h, m, eta, Q, HX, s)
        print("------------------------------double_effect_H2O_LiBr------------------------------")

        print_results_to_console(T, p, h, m, w, eta, Q, PP, s)

        if save_as_csv == "True":
            save_results_to_csv(T, p, h, m, w, eta, Q, PP, s)
    
    elif working_fluid == "NH3_H2O":
        T, p, h, m, w, eta, Q, PP, s = doubleEffect_model_NH3H2O.doubleEffect_model_NH3H2O(T, p, h, m, eta, Q, HX, s)
        print("------------------------------double_effect_NH3_H2O------------------------------")

        print_results_to_console(T, p, h, m, w, eta, Q, PP, s)

        if save_as_csv == "True":
            save_results_to_csv(T, p, h, m, w, eta, Q, PP, s)
    
    else:
        sys.exit("Working fluid is not defined!")

elif cycle_model == "DL":
    if working_fluid == "LiBr_H2O":
        T, p, h, m, w, eta, Q, PP, s = doubleLift_model_H2OLiBr.doubleLift_model_H2OLiBr(T, p, h, m, eta, Q, HX, s)
        print("------------------------------double_lift_H2O_LiBr------------------------------")

        print_results_to_console(T, p, h, m, w, eta, Q, PP, s)

        if save_as_csv == "True":
            save_results_to_csv(T, p, h, m, w, eta, Q, PP, s)
    
    elif working_fluid == "NH3_H2O":
        T, p, h, m, w, eta, Q, PP, s = doubleLift_model_NH3H2O.doubleLift_model_NH3H2O(T, p, h, m, eta, Q, HX, s)
        print("------------------------------double_lift_NH3_H2O------------------------------")

        print_results_to_console(T, p, h, m, w, eta, Q, PP, s)

        if save_as_csv == "True":
            save_results_to_csv(T, p, h, m, w, eta, Q, PP, s)
    
    else:
        sys.exit("Working fluid is not defined!")
