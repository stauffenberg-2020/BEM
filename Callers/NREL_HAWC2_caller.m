% Caller script
clearvars;clc;
N=3;
rho = 1.225;

op_pts = readmatrix('NREL5MWReferenceWindTurbinev50\data\operational_data.opt','FileType','text');
ae = readmatrix('NREL5MWReferenceWindTurbinev50\data\NREL_5MW_ae.txt');
pc = 'NREL5MWReferenceWindTurbinev50\data\NREL_5MW_pc.txt';

[profile_t_C, profile_AoA, profile_cL, profile_cD] = read_pc_file(pc);

wsp = op_pts(:,1);
pitch = op_pts(:,2);
rpm = op_pts(:,3);

sec_r = ae(2:size(ae,1),1);
sec_C = ae(2:size(ae,1),3);
sec_t_C = ae(2:size(ae,1),5);
sec_twist = [-13.3,-13.3,-13.3,-13.3,-11.48,-10.16,-9.011,-7.795,-6.544,-5.361,-4.188,-3.125,-2.319,-1.526,-0.863,-0.37,-0.106, 0]'; % to be read from *.htc file




induction_flag = 1; % 0 or 1, 0 = Induction Off, 1 = Induction On
tip_loss_flag = 1; % 0 or 1, 0 = Prandtl's tip loss correction Off, 1 = Prandtl's tip loss correction On
glauert_corr_flag = 2; % 0 or 1 or 2, 0 = Off, 1 = As per HANSEN eqn. 6.38, 2 = As per HANSEN eqn. 6.37 and as in Pyro

[output_details,output,BEM] = core_bem(rho, N, wsp, pitch, rpm, sec_r, sec_C, sec_t_C, sec_twist, profile_t_C, profile_AoA, profile_cL, profile_cD , induction_flag, tip_loss_flag, glauert_corr_flag);
% [output_details,output,BEM] = core_bem(rho, N, wsp(1), pitch(1), rpm(1), sec_r, sec_C, sec_t_C, sec_twist, profile_t_C, profile_AoA, profile_cL, profile_cD , induction_flag, tip_loss_flag, glauert_corr_flag);


%% Supporting functions

function [thickness, AoA, cL, cD] = read_pc_file(pc)
    fid = fopen(pc,'r');
    fid1 = fopen(pc,'r');
    
    nT = textscan(fid, '%f', 1, 'HeaderLines', 1, 'CollectOutput', 1);
    nT = nT{1,1}; % No. of tables
    entries = zeros(nT,1);
    
    j=1;
    for i=1:nT
        tline = textscan(fid, '%f%f%f', 1, 'HeaderLines', j, 'CollectOutput', 1);
        entries(i,1) = tline{1,1}(2);
        thickness(i,1) = tline{1,1}(3);
        if i==1
            section(i,1) = textscan(fid1, '%f%f%f%f%f%f', entries(i,1), 'HeaderLines', 3, 'CollectOutput', 1);
        else
            section(i,1) = textscan(fid1, '%f%f%f%f%f%f', entries(i,1), 'HeaderLines', 2, 'CollectOutput', 1);
        end
        j=entries(i,1)+1;
    end
    fclose(fid);
    fclose(fid1);
    % Linear interpolation to make sure all the data is in same size
    AoA = (-179:1:180)';
    for i=1:nT
        cL(:,i) = interp1(section{i,1}(:,1),section{i,1}(:,2),AoA);
        cD(:,i) = interp1(section{i,1}(:,1),section{i,1}(:,3),AoA);
    end
end

