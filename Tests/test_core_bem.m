folder = pwd;
filename = 'NREL_5MW.txt';   
file = fullfile(folder,filename);

parentDirectory = fileparts(cd);
addpath(genpath(parentDirectory));

[General, op_pts, Blade, ~] = read_turbine_file(file);

Blade.preflap = zeros(length(Blade.r),1);

[output_details, output, BEM] = core_bem(General, op_pts, Blade);

baseline_output = load(fullfile(folder,'baseline_output.mat'));

if round(output,4) == round(baseline_output.output,4)
    disp('Test Passed')
else
    throw(exception)
end
