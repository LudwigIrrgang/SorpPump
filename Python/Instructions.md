## Instructions Python-Tool SorpPump

# Needed Packages, Modules and Versions
The programm is developt within the Python 3.10.10 64-bit environment. In order to run the programm you need at least version 3.10 of Python.

The follwing librarys and packages have to be imported:
- time: Can be used for measuring the run time. Package is included but functionality has to be added later on.
- sys: Used for the sys.exit() function. (replacing error functionality from matlab).
- numpy: For the use of numpy arrays and other array related functions.
- pandas: Used for handeling data frames and exporting to .csv.
- math: extending the math functions of Python.
- scipy.optimize: Using different functions to find numerical solution of various problems.
- CoolProp.CoolProp: Importing the CoolProp library

# Basic operation

The model can be operated from the main.py file. To use the model two different approaches are provided. 
-   The user can choose to provide the input data within the main file. To do so the variable "external_input" in line 89 has to be "False". The data has to be provided in the lines 98 and following. The inputs varray based on the choosen model.
-   The other option is to provide the input data value by value with keyboard inputs. The required setting is "external_input = True" in line 89. The programm will than ask for the needed variables step by step.

The default setting is "external_input = True". If nothing is changed within the main.py file the cycle modelling can be performed without further knowlage of the python code. 

## More detailed description of the python-model

# main.py



