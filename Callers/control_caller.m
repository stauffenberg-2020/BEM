% Controller testing
clearvars;clc;close all;
file = '.\Data\NREL_5MW.txt';

[General, op_pts, Blade, CTRL] = read_turbine_file(file);
orig_pitch = op_pts.pitch;
Blade.preflap = zeros(length(Blade.r),1);

[~, output_orig, ~] = core_bem(General, op_pts, Blade);

[pitch, lambda, Cp] = read_Cp_Ct(CTRL.Cp_file);

MP = get(0,'MonitorPositions');
if size(MP,1) == 1
    SS = MP(1,:);
else
    SS = MP(2,:); % To have the results displayed better in bigger screen
end
figure('position',[SS(1)+SS(3)*0.05, SS(2)+SS(4)*0.2, SS(3)*0.9, SS(4)*0.6]);
subplot(1,2,1)
plotPitchVsLambda(pitch, lambda, Cp);
subplot(1,2,2)
plotCpVsLambda(pitch, lambda, Cp);

[op_pts.pitch, region] = core_pitch_ctr(op_pts.wsp, op_pts.rpm, General.rho, Blade.r(end), pitch, lambda, Cp, CTRL.Oin, CTRL.Orat, CTRL.Prat, 0.944);

[~, output, ~] = core_bem(General, op_pts, Blade);

%% Results Plotting
MP = get(0,'MonitorPositions');
if size(MP,1) == 1
    SS = MP(1,:);
else
    SS = MP(2,:); % To have the results better displayed in bigger screen
end
figure('position',[SS(1)+SS(3)*0.2, SS(2)+SS(4)*0.15, SS(4)*0.9, SS(4)*0.7])

plot(op_pts.wsp, orig_pitch,'Color','red','LineWidth',1.5,'LineStyle','--');
hold on
plot(op_pts.wsp, op_pts.pitch, 'Color','red','LineWidth',1.5,'LineStyle','-');
ax = gca;
ax.YColor = 'r';
ylabel('Pitch (deg)');

xl=xlim;
yl=ylim;
ids = find(logical(diff(region))==1);
x1 = [xl(1) op_pts.wsp(ids(1)) op_pts.wsp(ids(1)) xl(1)]; % wsp
[xl(1) op_pts.wsp(ids(1)) op_pts.wsp(ids(1)) xl(1)]; % wsp
y = [yl(1) yl(1) yl(2) yl(2)]; % pitch

x2 = [op_pts.wsp(ids(1)) op_pts.wsp(ids(2)) op_pts.wsp(ids(2)) op_pts.wsp(ids(1))]; % wsp
x3 = [op_pts.wsp(ids(2)) op_pts.wsp(ids(3)) op_pts.wsp(ids(3)) op_pts.wsp(ids(2))]; % wsp
x4 = [op_pts.wsp(ids(3)) xl(2) xl(2) op_pts.wsp(ids(3))]; % wsp

addaxis(op_pts.wsp, op_pts.rpm,'Color','blue','LineWidth',1.5,'LineStyle','--','DisplayName','RPM\_Orig');
hold on
addaxisplot(op_pts.wsp, op_pts.rpm,2,'Color','blue','LineWidth',1.5,'LineStyle','-','DisplayName','RPM\_Simple\_Ctr');
addaxislabel(2,'RPM');

addaxis(op_pts.wsp,output_orig(:,6)*0.944,'Color','black','LineWidth',1.5,'LineStyle','--','DisplayName','P\_Orig');
addaxisplot(op_pts.wsp,output(:,6)*0.944,3,'Color','black','LineStyle','-','LineWidth',1.5,'Color','k')
addaxislabel(3,'Power (kW)');


xlabel('Wind Speed (m/s)')
grid on

h(1) = patch(x1,y,'red','FaceAlpha',0.2,'DisplayName','const rpm');
h(2) = patch(x2,y,'blue','FaceAlpha',0.2,'DisplayName','varying rpm(max Cp)');
h(3) = patch(x3,y,'green','FaceAlpha',0.2,'DisplayName','transition');
h(4) = patch(x4,y,'yellow','FaceAlpha',0.2,'DisplayName','const rpm & power');

leg1 = legend('Pitch\_Orig','Pitch\_Simple\_Ctr','RPM\_Orig','RPM\_Simple\_Ctr','P\_Orig','P\_Simple\_Ctr','location','nw');
set(leg1,'FontSize',9);
ah1=axes('position',get(gca,'position'),'visible','off');
leg2 = legend(ah1, h(1:end),'location','se');set(leg2,'FontSize',9);
set(leg2,'FontSize',9);
title(leg2,'Regions')



