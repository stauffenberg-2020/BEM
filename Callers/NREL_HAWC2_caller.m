% Caller script
clearvars;clc;
N=3;
rho = 1.225;

op_data = readmatrix('G:\BEM\BEM\Data\NREL5MWRefTurb_v50\data\operational_data.opt','FileType','text');
ae = readmatrix('G:\BEM\BEM\Data\NREL5MWRefTurb_v50\data\NREL_5MW_ae.txt');
pc = 'G:\BEM\BEM\Data\NREL5MWRefTurb_v50\data\NREL_5MW_pc.txt';

op_pts.wsp = op_data(:,1);
op_pts.pitch = op_data(:,2);
op_pts.rpm = op_data(:,3);

sec_r = ae(2:size(ae,1),1);
sec_C = ae(2:size(ae,1),3);
sec_t_C = ae(2:size(ae,1),5);
sec_twist = [-13.3,-13.3,-13.3,-13.3,-11.48,-10.16,-9.011,-7.795,-6.544,-5.361,-4.188,-3.125,-2.319,-1.526,-0.863,-0.37,-0.106, 0]'; % to be read from *.htc file


BLD = organize_blade_data(sec_r, sec_C, sec_t_C, sec_twist, pc);

AeroFlags.induction = 1; % 0 or 1, 0 = Induction Off, 1 = Induction On
AeroFlags.tip_loss = 1; % 0 or 1, 0 = Prandtl's tip loss correction Off, 1 = Prandtl's tip loss correction On
AeroFlags.highCT = 2; % 0 or 1 or 2, 0 = Off, 1 = As per HANSEN eqn. 6.38, 2 = As per HANSEN eqn. 6.37

[output_details,output,BEM] = core_bem(rho, N, op_pts, BLD, AeroFlags);
% [output_details,output,BEM] = core_bem(rho, N, wsp(1), pitch(1), rpm(1), sec_r, sec_C, sec_t_C, sec_twist, profile_t_C, profile_AoA, profile_cL, profile_cD , induction_flag, tip_loss_flag, glauert_corr_flag);


%% Plotting & Comparisons

