#### Instructions Python-Tool SorpPump

### Needed Packages, Modules and Versions
The program is developed within the Python 3.10.10 64-bit environment. In order to run the program you need at least version 3.10 of Python.

The following libraries and packages have to be imported:
- time
- sys
- numpy
- pandas
- math
- scipy.optimize
- CoolProp.CoolProp: Importing the CoolProp library http://www.coolprop.org/
- ctREFPROP.ctREFPROP Implementation of the REFPROP Python Wrapper https://refprop-docs.readthedocs.io/en/latest/DLL/high_level.html

### Basic operation

The model can be operated from the main.py file. To use the model two different approaches are provided. 
-   The user can choose to provide the input data within the main file. To do so the variable "external_input" in line 86 has to be "False". The data has to be provided in the lines 94 and following. The inputs vary based on the chosen model.
-   The other option is to provide the input data value by value with keyboard inputs. The required setting is "external_input = True" in line 86. The program will than ask for the needed variables step by step.

The default setting is "external_input = True". If nothing is changed within the main.py file the cycle modelling can be performed without further knowledge of the python code. 

### More detailed description of the python-model

## main.py

After importing the required libraries and packages, the required modules "Fluids" and "Internal_Cycle_models" are imported. The containing functions are used to model the different fluids and cycle models.

The next step is the definition of two functions, which are used to handle the results. 
- "print_results_to_console(T, p, h, m, w, eta, Q, PP, s)" 
- "save_results_to_csv(T, p, h, m, w, eta, Q, PP, s)" 
Both functions use data frames to process the provided instances of the var() class. 
As the names imply the results are displayed in the console and/or saved in a .csv file.

Before the inputs are made, the var() class has to be defined. This class will be used to store the variables, similar to a matlab struct. The following types of variables are needed.

- temperature T [K]
- pressure p [Pa]
- enthalpy h [J/kg]
- mass flow m [kg/s]
- efficiency eta [-]
- Heat flow Q [W]
- Heat exchanger parameter HX [K]
- entropy and requirement s

# Necessary inputs


# Calculation



### Notation of the variables

The variables are named systematically. The first part, in front of the dot, describes the variable type. For example T.cond is a temperature value. 
After the variable type the fluid and location is defined. The description is a combination of the following options separated by "_".
    sol: solution of the two fluids
    ref: only the refrigerant
    abs: absorber
    des: desorber
    cond: condensator
    evap: evaporator
    pump: pump
    valve: expansion valve
    SHEX: Heat exchanger in solution line
    RHEX: Heat exchanger in refrigerant line
    in: in front of the component
    out: after the component

Is the description extended by "I" the variable is for a additional component or pressure level exceeding the base cycle.

    
