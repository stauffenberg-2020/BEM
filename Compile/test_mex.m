filename = 'NREL_5MW.txt';   

[General, op_pts, Blade, ~] = read_turbine_file(filename);

Blade.preflap = zeros(length(Blade.r),1);

[output_details, output, BEM] = core_bem(General, op_pts, Blade);

