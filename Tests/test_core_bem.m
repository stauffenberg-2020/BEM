file = '..\Data\NREL_5MW.txt';   
file = fullfile(pwd,file);

[General, op_pts, Blade, ~] = read_turbine_file(file);

Blade.preflap = zeros(length(Blade.r),1);

[output_details, output, BEM] = core_bem(General, op_pts, Blade);

if round(output(1,1),4) == 35.9070 % Entire output checking needs to be added
    disp('Test Passed')
else
    disp('Test Failed')
end
