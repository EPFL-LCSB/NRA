#### Matlab code for the Network Response Analysis (NRA) method required to replicate the results of:
#### "Constraint-based metabolic control analysis for rational strain engineering" by S. Tsouka, M. Ataman, T. Hameri, L. Miskovic, V. Hatzimanikatis.

This is a method to generate multiple metabolic engineering strategies given user defined options and limitations. Please refer to the documentation for basic instructions.



#### Requirements:

- matTFA, MATLAB toolbox for Thermodynamics-based Flux Balance Analyssis(TFA): avaliable at https://github.com/EPFL-LCSB/matTFA
- IBM CPLEX Studio, optimization suite, free for academic research
- Supplementary data files (users should place downloaded data in the "modeldata" subfolder): available at the Zenodo depository 10.5281/zenodo.4501189



#### Basic Instructions:

	1. The launch script for the method is “launch_script_NRA.m”. There you should adjust your local paths and options as you want and simply run it.

	2. The main NRA function is "run_NRA.m". Inside this function some important functions included are:
			i. The function "prepMCAcoefficientsForRBA.m" is the one that constructs the NRA matrices from the selected CCs. Most variable bounds and constraints (such as DG) are also added here.
			ii. The function "runConcViolationStudies.m" is the one that will perform the concentration violation study WITHOUT calculating alternatives for the violating metabolite sets(*).
			iii. The function "findAltControlSchemes.m" is the one that will calculate alternatives for strategies (sets of enzyme manipulations) for the given problem, up to the options.noAlt limit.


	(*)To generate alternatives for the violating metabolite species, run the script "Alternative_concs_for_violations.m". This script should be run INSIDE the generated folder from the concentration violation study as it uses its outputs.



#### Other remarks:

	1. The concentration violations are computed through the use of slack variables through the function "debugAdditionOfConstraints_W.m" (called inside of runConcViolationStudies.m).

	2. The alternative strategies are computed using integer cuts through the function "findDP_CONTROL_NEW.m" (called inside of findAltControlSchemes.m).



S.Tsouka et al. - Feb. 2021


