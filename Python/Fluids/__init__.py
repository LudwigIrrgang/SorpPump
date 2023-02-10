__all__ = ["LiBrSol", "NH3H2O", "H2O"]
import pandas as pd
from scipy.optimize import newton, fsolve
import math
import numpy
import CoolProp.CoolProp as CP
# global DFcal_H2O 
# DFcal_H2O = pd.read_csv("02_Test_module_structure\Fluids\H2O_cal.csv", sep = ";")
# global DFcal_LiBrSol
# DFcal_LiBrSol = pd.read_csv("02_Test_module_structure\Fluids\LiBrSol_cal.csv", sep=";")
