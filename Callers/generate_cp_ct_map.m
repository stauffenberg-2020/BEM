% Script for generating Cp and Ct map

% Caller script
clearvars;clc;
General.N=3;
General.rho = 1.225;

ae = readmatrix('G:\BEM\BEM\Data\NREL5MWRefTurb_v50\data\NREL_5MW_ae.txt');
pc = 'G:\BEM\BEM\Data\NREL5MWRefTurb_v50\data\NREL_5MW_pc.txt';

sec_r = ae(2:size(ae,1),1);
sec_C = ae(2:size(ae,1),3);
sec_t_C = ae(2:size(ae,1),5);
sec_twist = [-13.3,-13.3,-13.3,-13.3,-11.48,-10.16,-9.011,-7.795,-6.544,-5.361,-4.188,-3.125,-2.319,-1.526,-0.863,-0.37,-0.106, 0]'; % to be read from *.htc file

BLD = organize_blade_data(sec_r, sec_C, sec_t_C, sec_twist, pc);

General.induction = 1; % 0 or 1, 0 = Induction Off, 1 = Induction On
General.tip_loss = 1; % 0 or 1, 0 = Prandtl's tip loss correction Off, 1 = Prandtl's tip loss correction On
General.highCT = 2; % 0 or 1 or 2, 0 = Off, 1 = As per HANSEN eqn. 6.38, 2 = As per HANSEN eqn. 6.37

lambda = 5:0.5:10;
op_pts.pitch = (-5:0.5:14)';

op_pts.wsp = (3:2:18)';
map_pts.wsp = op_pts.wsp;
% for i=1:length(op_pts.pitch)
%     map_pts.pitch(1:length(op_pts.wsp),1) = op_pts.pitch(i);
%     for j=1:length(lambda)
%         map_pts.rpm = (lambda(j).*map_pts.wsp/sec_r(end)).*(30/pi);
%         
%         [~, output, ~] = core_bem(General, map_pts, BLD);
%         Cp(i,j) = max(output(:,7)); % The values are anyways same for a Lambda
%         Ct(i,j) = max(output(:,9));
%     end
% end

for i=1:length(op_pts.pitch)
    map_pts.pitch(1:length(op_pts.wsp),1) = op_pts.pitch(i);
    for j=1:length(lambda)
        map_pts.rpm = (lambda(j).*map_pts.wsp/sec_r(end)).*(30/pi);
        
        [~, output, ~] = core_bem(General, map_pts, BLD);
        Cp(i,j) = max(output(:,7)); % The values are anyways same for a Lambda
        Ct(i,j) = max(output(:,9));
    end
end

%% Plotting

close all;
figure()

contourf(op_pts.pitch,lambda,Cp','ShowText','on')
xlabel('Pitch (deg)');
ylabel('Lambda (TSR)');
hcb = colorbar;
hcb.Title.String = 'Cp';

%%
figure()
for i=1:4:length(op_pts.pitch)
    plot(lambda,Cp(:,i),'DisplayName',num2str(op_pts.pitch(i)))
%     ax = plotcontour(lambda,Cp(:,i),op_pts.pitch(i));
    hold on
    [max_Cp,idx] = max(Cp(:,i));
    plot(lambda(idx),max_Cp,'r*','HandleVisibility','off');
end
grid on
leg=legend('location','nw');
title(leg,'Pitch (deg)');
ylim([0 0.593])
ylabel('Cp')
xlabel('Lambda (TSR)')

%%
figure()
plotCpVsLambda(lambda,op_pts.pitch(1:6:39)',(Cp(:,1:6:39)));

function [CM,cc] = plotCpVsLambda(inp_x,inp_z,inp_y)
    
    % inp_x: 1xM-vector with x-coordinates
    % inp_z: 1xN-vector with z-coordinates (the 'text-Info')
    % inp_y: MxN-matrix with datapoints

    [X,Y]=meshgrid(inp_x,inp_z);
    [X_inp,Y_inp]=meshgrid(inp_x,inp_z);
    Z=interp2(X_inp,Y_inp,inp_y',X,Y);
    [CM, cc] = contour(X,Z,Y,inp_z,'ShowText','on');
    ylim([0 0.593]);
    grid on
    xlabel('Lambda (TSR)')
    ylabel('Cp')
    hcb = colorbar;
    hcb.Title.String = 'Pitch (deg)';
end