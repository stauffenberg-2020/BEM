files = dir(pwd);
dirFlags = [files.isdir];
subFolders = files(dirFlags);
subFolderNames = {subFolders(3:end).name}; 
if isempty(subFolderNames)
    cd ..
elseif ~any(strcmp(subFolderNames,'Tests'))
    cd ..
end
addpath(genpath(pwd))

cd Compile/
warning('off','all');

licen = evalc("license checkout MATLAB_Coder");
licens = str2double(extractAfter(licen,'='));

if licens == 1
    % Build library
    codegen -report -config:lib bem_calc.m -args {coder.typeof('NREL_5MW.txt',[1,inf])}
    disp('Library build successful')

    % Build MEX
    clear functions
    codegen -report bem_calc.m -args {'NREL_5MW.txt'} -test test_core_bem
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
else
    disp('Coder license not available')
end
