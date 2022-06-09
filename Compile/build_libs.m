codegen -report -config:lib bem_calc.m -args {'G:\BEM\BEM\Data\NREL_5MW.txt'}

clear functions
codegen -report bem_calc.m -args {'G:\BEM\BEM\Data\NREL_5MW.txt'} -test LPE_example

