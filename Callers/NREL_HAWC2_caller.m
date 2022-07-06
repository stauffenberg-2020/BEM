% Caller script
clearvars;clc;

filename = 'NREL_5MW.txt';
[General, op_pts, BLD, ~] = read_turbine_file(filename);

General.highCT = 1; % 0 or 1 or 2, 0 = Off, 1 = As per HANSEN eqn. 6.38, 2 = As per HANSEN eqn. 6.37
BLD.preflap = zeros(length(BLD.r),1);

plot(op_pts.wsp,op_pts.pitch)
hold on
op_pts.pitch(9:end) = [3.83,6.60,8.7,10.45,12.06,13.54,14.92,16.23,17.47,18.70,19.94,21.18,22.35,23.47]'; % Values from https://www.nrel.gov/docs/fy09osti/38060.pdf
plot(op_pts.wsp,op_pts.pitch)

[output_details, output, BEM] = core_bem(General, op_pts, BLD);

op_data = readmatrix('/Data/NREL5MWRefTurb_v50/data/operational_data.opt','FileType','text');

%% Plotting & Comparisons
% close all
figure()
plot(op_pts.wsp,output(:,6))
hold on
plot(op_pts.wsp,op_data(:,4))
grid on
xlabel('Wind Speed (m/s)')
ylabel('Aerodynamic Power (kW)')
legend('BEM','NREL-H2','location','nw')

figure()
plot(op_pts.wsp,output(:,8))
hold on
plot(op_pts.wsp,op_data(:,5))
grid on
xlabel('Wind Speed (m/s)')
ylabel('Aerodynamic Thrust (kN)')
legend('BEM','NREL-H2','location','nw')


