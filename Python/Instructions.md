#### Instructions Python-Tool SorpPump

### Needed Packages, Modules and Versions
The program is developed within the Python 3.10.10 64-bit environment. In order to run the program you need at least version 3.10 of Python.

The following libraries and packages have to be imported:
- time: Can be used for measuring the run time. Package is included but functionality has to be added later on.
- sys: Used for the sys.exit() function. (replacing error functionality from matlab).
- numpy: For the use of numpy arrays and other array related functions.
- pandas: Used for handling data frames and exporting to .csv.
- math: extending the math functions of Python.
- scipy.optimize: Using different functions to find numerical solution of various problems.
- CoolProp.CoolProp: Importing the CoolProp library
- ctREFPROP.ctREFPROP


Implementation of the REFPROP Python Wrapper
https://refprop-docs.readthedocs.io/en/latest/DLL/high_level.html
And definition of the used Unit set


### Basic operation

The model can be operated from the main.py file. To use the model two different approaches are provided. 
-   The user can choose to provide the input data within the main file. To do so the variable "external_input" in line 89 has to be "False". The data has to be provided in the lines 98 and following. The inputs vary based on the chosen model.
-   The other option is to provide the input data value by value with keyboard inputs. The required setting is "external_input = True" in line 89. The program will than ask for the needed variables step by step.

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

# Data input

The input data can be provided in two different ways. Directly in the main.py file are step by step as keyboard inputs in the Terminal. 

To use the input provided in the main.py file the variable "external_input" in line 92 has to be set to the Boolean False. If it is set to True the program will ask for each input individually.

If external_input is set to False the input data provided in the lines 100 to 152 are used for calculation. The possible options and short descriptions for each parameter are provided within the code. The Variable names and the syntax of the inputs is the same as the ones for external inputs.

Is external_input set to True the program will request each input. 
The first two parameters define which fluid and cycle model shall be modeled. 
- "working_fluid" defines which mixture is used. It can be set to two different strings: NH3_H2O and LiBr_H2O. All inputs are case sensitive. 
- "cycle_model" chooses the desired cycle configuration. The possible options are the strings: base, DEK and DLK. DEK = Double Effect Cycle and DLK = Double Lift Cycle.

With the parameter "save_as_csv" the user decides if the results will be saved as .csv or will only be printed in the Terminal. The options are "True" and "False". The input is handled as string not as Boolean. 

After the general configurations the values for the input data are defined. 
At first the information needed for all cycle types and working fluids are defined. 
The following variables have to be assigned with a value:
    T.evap = Temperature which is needed in the evaporator
    T.sol_abs_out = Temperature of the solution leaving the absorber (lowest pressure level)
    T.sol_des_out = Temperature of the solution leaving the desorber (one pressure level higher than T.sol_abs_out)
    T.cond = The condensation temperature.
    T.ext_cond_in = The starting temperature of the external cooling solution

    s.requirement = Which heat flow will be provided? (Q_evap: heat flow at the evaporator or Q_des: heat flow at the desorber)

    Q.dec = The heat flow value based on s.requirement

    eta.pump = The efficiency of the pump

    HX.T_PP_SHEX = Heat exchanger parameter of the solution heat exchanger
    HX.T_PP_RHEX = Heat exchanger parameter of the refrigerant heat exchanger
    HX.T_PP_cond = Heat exchanger parameter of the condensator

Is the cycle_model set to DEK further inputs are needed:
    T.sol_des_outI = Temperature of the solution leaving the desorber (highest pressure level)

    HX.T_PP_SHEXI = Heat exchanger parameter of the solution heat exchanger in the high pressure cycle
    HX.T_PP_cond_int = Heat exchanger parameter of the second condenser

Is the cycle_model set to DLK further inputs are needed:
    T.sol_abs_outI = Temperature of the solution leaving the absorber (middle pressure level)
    T.sol_des_outI = Temperature of the solution leaving the desorber (highest pressure level)

    HX.T_PP_SHEXI = Heat exchanger parameter of the solution heat exchanger in the high pressure cycle

For the mixture of ammonia and water further inputs are required:
    All there cycle types require the following parameter:
        HX.dT_ref_des = Temperature difference of the refrigerant in the desorber
    The DLK and DEK require the following additional parameter
        HX.dT_ref_desI = Temperature difference of the refrigerant in the desorber of the highest pressure level

# Calculation of the model

The defined "cycle_model" and "working_fluid" determine which function will be called. 
All functions have the same input parameters:
    T, p, h, m, eta, Q, HX, s
The parameter are of the class-type as described earlier.

Within the functions the calculations will take place. After the calculations have finished the results are provided in the same classes as the inputs, extended by the calculated values. The returned parameters are:
    T, p, h, m, w, eta, Q, PP, s
With "w" the mass fraction of the solutions and "PP" the process parameter or post processing data.

The results are printed to the terminal and can be saved as .csv file. Additional functionality, like plotting graphs, can be added at this state.

### Notation of the variables

The variables are named systematical. The first part, in front of the dot, describes the variable type. For example T.cond is a temperature value. 
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

    
