% Turbine caller
clearvars;clc;
file = 'G:\BEM\BEM\Data\NREL_5MW.txt';
[General, op_pts, Blade] = read_turbine_file(file);

Blade.preflap = zeros(length(Blade.r),1);

AeroFlags.induction = General.induction; % 0 or 1, 0 = Induction Off, 1 = Induction On
AeroFlags.tip_loss = General.tip_loss; % 0 or 1, 0 = Prandtl's tip loss correction Off, 1 = Prandtl's tip loss correction On
AeroFlags.highCT = General.highCT; % 0 or 1 or 2, 0 = Off, 1 = As per HANSEN eqn. 6.38, 2 = As per HANSEN eqn. 6.37

[output_details,output,BEM] = core_bem(General, op_pts, Blade);