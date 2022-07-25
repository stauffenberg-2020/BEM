% Independent Blade calc test
clearvars;clc;close all;

filename = 'NREL_5MW.txt';
[General, op_pts, BLD, ~] = read_turbine_file(filename);

General.induction = 1;
General.highCT = 1; % 0 or 1 or 2, 0 = Off, 1 = As per HANSEN eqn. 6.38, 2 = As per HANSEN eqn. 6.37
BLD.preflap = zeros(length(BLD.r),1);
op_pts.shear(1) = 0.3;

plot(op_pts.wsp,op_pts.pitch)
hold on
op_pts.pitch(9:end) = [3.83,6.60,8.7,10.45,12.06,13.54,14.92,16.23,17.47,18.70,19.94,21.18,22.35,23.47]'; % Values from https://www.nrel.gov/docs/fy09osti/38060.pdf
plot(op_pts.wsp,op_pts.pitch)

tic
[output_details, output, BEM] = core_bem_orig(General, op_pts, BLD);
toc
tic
[~, output_ind, BEM_ind] = core_bem(General, op_pts, BLD);
toc

op_data = readmatrix('/Data/NREL5MWRefTurb_v50/data/operational_data.opt','FileType','text');

%% Plotting & Comparisons
% close all
figure()
plot(op_pts.wsp,output(:,6))
hold on
plot(op_pts.wsp,output_ind(:,6),'--')
plot(op_pts.wsp,op_data(:,4))
grid on
xlabel('Wind Speed (m/s)')
ylabel('Aerodynamic Power (kW)')
legend('Orig BEM No Ind','New BEM No Ind','NREL-H2','location','nw')

figure()
plot(op_pts.wsp,output(:,8))
hold on
plot(op_pts.wsp,output_ind(:,8),'--')
plot(op_pts.wsp,op_data(:,5))
grid on
xlabel('Wind Speed (m/s)')
ylabel('Aerodynamic Thrust (kN)')
legend('Orig BEM No Ind','New BEM No Ind','NREL-H2','location','nw')

%% Plotting shear profile
HH = 90;
wsp_id = 1;
MP = get(0,'MonitorPositions');
if size(MP,1) == 1
    SS = MP(1,:);
else
    SS = MP(2,:); % To have the results better displayed in bigger screen
end
figure('position',[SS(1)+SS(3)*0.1, SS(2)+SS(4)*0.15, SS(4)*1.4, SS(4)*0.7]);
subplot(1,2,1)
heights = [flip(HH-BEM.r); HH; HH+BEM.r];
sheared_wsp = op_pts.wsp(wsp_id).*(heights./HH).^(op_pts.shear(wsp_id));
surf(heights, sheared_wsp, repmat(heights,1,length(sheared_wsp)),'EdgeColor','flat','FaceColor','interp','LineWidth',2)
view(90,0)
hold on
plot3([HH HH],[op_pts.wsp(wsp_id) op_pts.wsp(wsp_id)],[min(heights) max(heights)],'LineWidth',1.5,'color','k')
xline(HH,'displayname','HH')
ylabel('Wind Speed variation');
zlabel('Rotor height');
% hcb = colorbar;
% title(hcb,'Wind Shear');

subplot(1,2,2)
fig = plot_over_disc(BEM_ind, wsp_id, 'wsp_sh');
