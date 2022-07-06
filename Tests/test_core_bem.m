filename = 'NREL_5MW.txt';
cd ..
addpath(genpath(pwd))

[General, op_pts, Blade, ~] = read_turbine_file(filename);

Blade.preflap = zeros(length(Blade.r),1);

[output_details, output, BEM] = core_bem(General, op_pts, Blade);

try 
    load('NREL_5MW_output.mat');
catch
    disp('Unable to load baseline data')
end

if max(abs((output./baseline_output)-1)*100,[],'all')<0.01
    disp('Test Passed, Changes are within 0.01% w.r.t baseline')
else
    disp('Results variation are more than 0.01%')
    throw(exception)
end
