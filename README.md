# matricesDeContactoMDP

The Python codes presented here have been used to compute the results presented in the following scientific article: https://www.researchgate.net/publication/346701109_An_epidemic_model_for_COVID-19_transmission_in_Argentina_Exploration_of_the_alternating_quarantine_and_massive_testing_strategies

The authors share this software "as is", for reference and reproducibility purposes. Any query can be directed to the corresponding author.

About the code:

deterministicEqs_v5.2.py: This program calculates variables of interest as a function of time of an extended SEIRD epidemic model. The system of equations is solved numerically using the explicit Euler method. Many of the model parameters are read from the attached folders, including population contact matrices (read below). Inside ParamEtar folder, the initial value of each epidemic compartment (also divided by age cohorts) is set in "input.txt". In addition, a list with different parameter values can be written to the sim_list.txt file. The program calculates the results corresponding to each item in this queue and saves them in different folders. In addition, an extra folder called "res_meta_det" will be created where final values of some variables of interest are saved as a function of the parameters provided in the list.  

stochasticEqs_v5.1.py: This program solves the same equations as the previous one but in its stochastic version. It does not have implemented the reading of the sim_list.txt file yet. The parameters provided there must be entered manually in the body of the program. In turn, it must be indicated how many realizations will be made for this set of values, through the variable "nrea" (line 14). For optimization purposes, some of its routines have been written in Cython (C-Extensions for Python). Therefore, follow the steps below before running it:

1 - Install Cython. Documentation: https://cython.readthedocs.io/en/latest/

2- From the console, compile setup.py:
python setup.py build_ext --inplace

3- Run stochasticEqs_v5.1.py

About the contact matrices:

Contact matrices correspond to the inhabitants of the city of Mar del Plata, Argentina. They were estimated from official census data (1) and the inference method developed by L Fumanelli et.al. (2). 21 age groups of ages 0-4, 5-9, ..., 100-104 are used. All matrices are normalized.

(1) Minnesota Population Center. Integrated Public Use Microdata Series, International: Version 7.2 [dataset]; Minneapolis, MN: IPUMS, 2019.  https://doi.org/10.18128/D020.V7.2

(2) Fumanelli, L., Ajelli, M., Manfredi, P., Vespignani, A., & Merler, S. (2012). Inferring the structure of social contacts from demographic data in the analysis of infectious diseases spread. PLoS Comput Biol, 8(9), e1002673. https://doi.org/10.1371/journal.pcbi.1002673
