folder = pwd;
filename = 'NREL_5MW.txt';   
file = fullfile(folder,filename);

parentDirectory = fileparts(cd);
addpath(genpath(parentDirectory));

[General, op_pts, Blade, ~] = read_turbine_file(file);

Blade.preflap = zeros(length(Blade.r),1);

[output_details, output, BEM] = core_bem(General, op_pts, Blade);

baseline_output = load('baseline_output');

if output == baseline_output.output % Entire output checking needs to be added
    disp('Test Passed')
else
    throw(exception)
end
