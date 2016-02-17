Instructions for figures from the Theisen CBB paper.

MATLAB v.7.12.0

Fig 1) Freehand figure in PowerPoint

Fig 2) 
-Execute Modeler_noPPP.m
	+ First Enter a suitable name in the quotes on line 132, for example: 'noPPP_WPi_300', to indicate no PPP with phosphate, 300 models
- To model without phosphate, comment out lines 61-68 and enter a suitable name like 'noPPP_NoPi_300'
- Insert the two model result names into the quotes in curly brackets on line 7 of Plot_Robustness_multi.m and execute
	+ brackets should look like: {'noPPP_WPi_300'; 'noPPP_NoPi_300'} or equivalent

Fig 3) 
-To model TPI with regulation, uncomment lines 128 & 129 & ensure that lines 61-68 are uncommented to include phosphate in the model. Choose a suitable name for line 132 such as 'noPPP_WPi_WTpiReg_300'.
- Insert the two model result names into the quotes in curly brackets on line 7 of Plot_Robustness_multi.m and execute
	+ Change line 34 of Plot_Robustness_multi.m to: for Count = 4; 
	+ Uncomment line 41.
	+ brackets should look like: {'noPPP_WPi_300'; 'noPPP_WPi_WTpiReg_300'} or equivalent

Fig 4)
-To model TPI with 10% G6P shunt, run Modeler_oPPP 
	+ First Enter a suitable name in the quotes on line 132, for example: '10_pct_oPPP_300', to indicate no PPP with phosphate, 300 models
-To model TPI with 30% G6P shunt:
	+ Change the value of the integer on line 96 to 3, change the name on line 132 to for example '30_pct_oPPP_300'
- Insert the three model result names into the quotes in curly brackets on line 7 of Plot_Robustness_multi.m and execute
	+ Revert all changes from Fig 3 (reload from .zip)
	+ brackets should look like: {'noPPP_WPi_300'; '10_pct_oPPP_300'; '30_pct_oPPP_300'} or equivalent
	+ uncomment line 44, make other changes so that plot function extends to line 44 and includes all 3 model plots.

Fig 5)
- Run FindFluxesOne.m with the appropriate appropriate model results with 'noPPP_WPi_300' & '10_pct_oPPP_300' & '30_pct_oPPP_300' in the quotes on line 6.
- Variable names:
	+ A: Average final RuBisco (inital is 10) flux after perturbation of 12 enzymes in CBB
	+ B: standard dev of A
	+ C: Average final G6PDH flux (inital is 1 or 3) after perturbation of 12 enzymes in CBB (may not exist if model has no G6PDH reaction e.g. noPPP)
	+ D: standard dev of C

- Calculate net CO2 fixation rate as A-C (or A if C does not exist), percent change will be based on difference from 10, 9 or 7 as appropriate

