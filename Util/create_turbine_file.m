% Script for creating turbine file from HAWC2 setup
clearvars;close all;clc;

filename = 'Data/NREL_5MW.txt';

General.N=3;
General.rho = 1.225;
General.induction = 1; % 0 or 1, 0 = Induction Off, 1 = Induction On
General.tip_loss = 1; % 0 or 1, 0 = Prandtl's tip loss correction Off, 1 = Prandtl's tip loss correction On
General.highCT = 2; % 0 or 1 or 2, 0 = Off, 1 = As per HANSEN eqn. 6.38, 2 = As per HANSEN eqn. 6.37

Control.CTR_flag = 0; % Simple pitch controller flag (0/1). If 1, pitch from OPERATIONAL_SET_POINTS will not be used and will be calculated
Control.Cp_file = 'Data/Cp_data.txt'; % Path to Cp matrix
Control.Ct_file = 'Data/Ct_data.txt'; % Path to Ct matrix
Control.Prat = 5000; % Rated Power in kW
Control.Vin = 3; % Cut in wind speed (m/s)
Control.Vrat = 11.4; % Rated wind speed (m/s)
Control.Vout = 25; % Cut out wind speed (m/s)
Control.Oin = 6.9; % Cut in rotor rpm
Control.Orat = 12.1; % Rated rotor rpm

op_data = readmatrix('./Data/NREL5MWRefTurb_v50/data/operational_data.opt','FileType','text');
ae = './Data/NREL5MWRefTurb_v50/data/NREL_5MW_ae.txt'; % ae data file
htc = './Data/NREL5MWRefTurb_v50/htc/NREL_5MW_reference_wind_turbine.htc'; % htc file
pc = 'Data/NREL5MWRefTurb_v50/data/NREL_5MW_pc.txt'; % Aerofoil data file

[BLD.r, BLD.C, BLD.t_C, nsec] = read_ae_file(ae); 
[BLD.AeroTwist] = read_aero_twist(htc, nsec);

Data.GENERAL = General;
Data.OPERATIONAL_SET_POINTS = op_data;
Data.BLADE_DETAILS = BLD;
Data.AEROFOIL_FILE = pc;
Data.CONTROL = Control;

write_turbine_file(Data,filename);

%% Supporting funtions

function write_turbine_file(Data,filename)
    fileID = fopen(filename,'w');
    header_titles = fieldnames(Data);
    
    for i=1:length(header_titles)
        if strncmpi(header_titles{i}, 'GENERAL',7)
            fn = fieldnames(Data.(header_titles{i}));
            fprintf(fileID,'# %s\n',header_titles{i});
            fprintf(fileID,'%d; %% No. of blades\n',Data.(header_titles{i}).(fn{1}));
            fprintf(fileID,'%0.3f; %% Density of air\n',Data.(header_titles{i}).(fn{2}));
            fprintf(fileID,'%d; %% Induction off/on = 0/1\n',Data.(header_titles{i}).(fn{3}));
            fprintf(fileID,'%d; %% Prandtl''s tip loss correction off/on = 0/1\n',Data.(header_titles{i}).(fn{4}));
            fprintf(fileID,'%d; %% Glauert''s high CT correction off/per HANSEN eqn 6.38/per HANSEN eqn 6.37 = 0/1/2\n',Data.(header_titles{i}).(fn{5}));
            fprintf(fileID,'\n');
        elseif strncmpi(header_titles{i}, 'OPERATIONAL_SET_POINTS',22)
            fprintf(fileID,'# %s\n',header_titles{i});
            fprintf(fileID,'wsp\tshear\tpitch\trpm\n');
            for j=1:size(Data.(header_titles{i}),1)
                fprintf(fileID,'%d\t%0.1f\t%0.3f\t%0.2f\n',Data.(header_titles{i})(j,1),0,Data.(header_titles{i})(j,2),Data.(header_titles{i})(j,3));
            end
            fprintf(fileID,'\n');
        elseif strncmpi(header_titles{i}, 'BLADE_DETAILS',13)
            fprintf(fileID,'# %s\n',header_titles{i});
            fprintf(fileID,'r\tC\tt_C\tAeroTwist\n');
            for j=1:length(Data.(header_titles{i}).r)
                fprintf(fileID,'%0.3f\t%0.3f\t%d\t%0.3f\n',Data.(header_titles{i}).r(j,1),Data.(header_titles{i}).C(j,1),Data.(header_titles{i}).t_C(j,1),Data.(header_titles{i}).AeroTwist(j,1));
            end
            fprintf(fileID,'\n');
        elseif strncmpi(header_titles{i}, 'AEROFOIL_FILE',13)
            fprintf(fileID,'# %s\n',header_titles{i});
            fprintf(fileID,'%s\n',Data.(header_titles{i}));
            fprintf(fileID,'\n');
        elseif strncmpi(header_titles{i}, 'CONTROL',7)
            fn = fieldnames(Data.(header_titles{i}));
            fprintf(fileID,'# %s\n',header_titles{i});
            fprintf(fileID,'%d; %% Simple pitch controller flag (0/1). If 1, pitch from OPERATIONAL_SET_POINTS will not be used and instead calculated\n',Data.(header_titles{i}).(fn{1}));
            fprintf(fileID,'%s; %% Path to Cp matrix\n',Data.(header_titles{i}).(fn{2}));
            fprintf(fileID,'%s; %% Path to Cp matrix\n',Data.(header_titles{i}).(fn{3}));
            fprintf(fileID,'%d; %% Rated Power in kW\n',Data.(header_titles{i}).(fn{4}));
            fprintf(fileID,'%0.1f; %% Cut in wind speed (m/s)\n',Data.(header_titles{i}).(fn{5}));
            fprintf(fileID,'%0.1f; %% Rated wind speed (m/s)\n',Data.(header_titles{i}).(fn{6}));
            fprintf(fileID,'%0.1f; %% Cut out wind speed (m/s)\n',Data.(header_titles{i}).(fn{7}));
            fprintf(fileID,'%0.1f; %% Cut in rotor rpm\n',Data.(header_titles{i}).(fn{8}));
            fprintf(fileID,'%0.1f; %% Rated rotor rpm\n',Data.(header_titles{i}).(fn{9}));
        end
        
    end
       
    fclose(fileID);
end
