filename = 'NREL_5MW.txt';
cd ..
addpath(genpath(pwd))

[General, op_pts, Blade, ~] = read_turbine_file(filename);

Blade.preflap = zeros(length(Blade.r),1);

[output_details, output, BEM] = core_bem(General, op_pts, Blade);

try 
    baseline_output = load('baseline_output.mat');
catch
    disp('Unable to load baseline data')
end

if round(output,4) == round(baseline_output.output,4)
    disp('Results matching upto 4th decimal, Test Passed')
else
    disp('Results not matching')
    throw(exception)
end
