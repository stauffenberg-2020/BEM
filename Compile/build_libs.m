% Build library
cd Compile\
warning('off','all');
codegen -report -config:lib bem_calc.m -args {'NREL_5MW.txt'}
disp('Library build successful')

% Build MEX
clear functions
codegen -report bem_calc.m -args {'NREL_5MW.txt'} -test test_mex
disp('MEX build successful')

% Build EXE
cfg = coder.config( "exe", "ecoder", false );
cfg.CustomInclude = "C_Files";
cfg.CustomSource = "C_Files/main.c";
cfg.EnableAutoCommit = false;
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;
codegen -config cfg bem_calc.m -args {coder.typeof('NREL_5MW.txt',[1,inf])}
disp('EXE build successful')
warning('on','all')

