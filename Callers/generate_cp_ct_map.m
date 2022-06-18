% Script for generating Cp (and Ct) map

clearvars;close all;clc;

ae = '.\Data\NREL5MWRefTurb_v50\data\NREL_5MW_ae.txt'; % ae data file
htc = '.\Data\NREL5MWRefTurb_v50\htc\NREL_5MW_reference_wind_turbine.htc'; % htc file
pc = '.\Data\NREL5MWRefTurb_v50\data\NREL_5MW_pc.txt'; % Aerofoil data file

[BLD.r, BLD.C, BLD.t_C, nsec] = read_ae_file(ae);
[BLD.AeroTwist] = read_aero_twist(htc, nsec);
[BLD.pro_t_C, BLD.pro_AoA, BLD.pro_cL, BLD.pro_cD] = read_pc_file(pc);
BLD.preflap = zeros(length(BLD.r),1);

General.N=3;
General.rho = 1.225;
General.induction = 1; % 0 or 1, 0 = Induction Off, 1 = Induction On
General.tip_loss = 1; % 0 or 1, 0 = Prandtl's tip loss correction Off, 1 = Prandtl's tip loss correction On
General.highCT = 2; % 0 or 1 or 2, 0 = Off, 1 = As per HANSEN eqn. 6.38, 2 = As per HANSEN eqn. 6.37

lambda = unique([2:0.2:3.8 4.25:0.25:6.75 7.5:0.5:8 8:1:20]); % Lambda (TSR) range, 2:1:20
pitch = unique([-5:0.25:5 5:0.5:10 10:1:25])'; % Pitch range, -5:0.5:25

wsp = 8; % Just a temp variable for deriving rpm corresponding to lambda

Cp = zeros(length(pitch),length(lambda));
Ct = zeros(length(pitch),length(lambda));
for i=1:length(pitch)
    op_pts.wsp(1:length(lambda),1) = wsp;
    op_pts.pitch(1:length(lambda),1) = pitch(i);
    op_pts.rpm = (lambda'*wsp/BLD.r(end)).*(30/pi);
    [~, output, ~] = core_bem(General, op_pts, BLD);
    Cp(i,:) = output(:,7);
    Ct(i,:) = output(:,9);
end
%
path = '.\Data\';
write_Cp_Ct(lambda,pitch,Cp,Ct,path);
%% Plotting Cp & Ct

MP = get(0,'MonitorPositions');
if size(MP,1) == 1
    SS = MP(1,:);
else
    SS = MP(2,:); % To have the results better displayed in bigger screen
end
figure('position',[SS(1)+SS(3)*0.05, SS(2)+SS(4)*0.2, SS(3)*0.9, SS(4)*0.6]);
subplot(1,2,1)
[C,h] = plotPitchVsLambda(pitch, lambda, Cp);
subplot(1,2,2)
[CM,cc] = plotCpVsLambda(pitch, lambda, Cp);

%% Upscaling the CP and CT data for better plotting purposes
[lambda_up, pitch_up, Cp_up, Ct_up] = upscale(lambda, pitch, Cp, Ct, 200, 500);

positive_Cp_data = max(min(Cp_up, 0.593), 0);
flag = ceil(positive_Cp_data);
Cp_up = Cp_up.*flag;
Ct_up = Ct_up.*flag;

%% Cp and Ct plot
MP = get(0,'MonitorPositions');
if size(MP,1) == 1
    SS = MP(1,:);
else
    SS = MP(2,:); % To have the results better displayed in bigger screen
end
figure('position',[SS(1)+SS(3)*0.05, SS(2)+SS(4)*0.2, SS(3)*0.9, SS(4)*0.6])
subplot(1,2,1)
levels = unique([flip(min(Cp_up(:)):0.001:0.01) 0.01:0.05:max(Cp_up(:))]); % Manually selecting levels for better plotting
[C2,h2] = contourf(lambda_up', pitch_up, Cp_up, levels,'LineStyle','--');
hold on
xlabel('Lambda (TSR)');
ylabel('Pitch (deg)');
title('Cp map');

Lvls2 = h2.LevelList;
h2.LevelList = Lvls2(Lvls2 > 0); 
h2.LevelList = round(h2.LevelList,3);
v2 = unique([Lvls2(1) Lvls2(11:end)]); % Manually selecting levels for labeling
clabel(C2,h2,v2);

[x,y]=find(ismember(Cp_up,max(Cp_up(:))));
plot(lambda_up(y),pitch_up(x),'r-p')
txt = sprintf('  Cp_m_a_x = %0.3f',max(Cp_up(:)));
text(lambda_up(y),pitch_up(x),txt,'Color','red');
line([lambda_up(1),lambda_up(y)],[pitch_up(x),pitch_up(x)],'Color','red','LineStyle','--')
line([lambda_up(y),lambda_up(y)],[pitch_up(1),pitch_up(x)],'Color','red','LineStyle','--')
hcb2 = colorbar;
hcb2.Title.String = 'Cp';

for j=1:length(lambda_up)
    if max(Cp_up(:,j)) == 0
        I_up(j) = I_up(j-1);
    else
        [~,I_up(j)] = max(Cp_up(:,j));
    end
end
plot(lambda_up,pitch_up(I_up),'-r','LineWidth',0.5);


subplot(1,2,2)
levels = unique([flip(min(Ct_up(:)):0.02:0.1) 0.2:0.2:max(Ct_up(:))]); % Manually selecting levels for better plotting
[C3,h3] = contourf(lambda_up', pitch_up, Ct_up, levels,'LineStyle','--');
hold on
Lvls3 = h3.LevelList;
h3.LevelList = Lvls3(Lvls3 > 0); 
h3.LevelList = round(h3.LevelList,2);
v3 = unique([Lvls3(1) Lvls3(6:end)]); % Manually selecting levels for labeling
clabel(C3,h3,v3);
hcb3 = colorbar;
hcb3.Title.String = 'Ct';
xlabel('Lambda (TSR)');
ylabel('Pitch (deg)');
title('Ct map (corresponding to Cp map)');


%% Supporting functions 

function write_Cp_Ct(lambda,pitch,Cp,Ct,path)
    % Writing the outputs into a text file
    cp_file = plus(path,"Cp_data.txt");
    ct_file = plus(path,"Ct_data.txt");
    fileID1 = fopen(cp_file,'w');
    fileID2 = fopen(ct_file,'w');
    
    fprintf(fileID1,'%d \t %d\n',length(pitch),length(lambda));
    fprintf(fileID2,'%d \t %d\n',length(pitch),length(lambda));
    
    for i=1:length(pitch)+1
        if i==1
            for j=1:length(lambda)+1
                if j==1
                    fprintf(fileID1,'            \t');
                    fprintf(fileID2,'            \t');
                else
                    fprintf(fileID1,'%8.4f \t',lambda(j-1));
                    fprintf(fileID2,'%8.4f \t',lambda(j-1));
                end
                if j==(length(lambda)+1)
                    fprintf(fileID1,'\n');
                    fprintf(fileID2,'\n');
                end
            end
        else
            for j=1:length(lambda)+1
                if j==1
                    fprintf(fileID1,'%8.3f \t',pitch(i-1));
                    fprintf(fileID2,'%8.3f \t',pitch(i-1));
                else
                    fprintf(fileID1,'%8.4f \t',Cp(i-1,j-1));
                    fprintf(fileID2,'%8.4f \t',Ct(i-1,j-1));
                end
                if j==(length(lambda)+1)
                    fprintf(fileID1,'\n');
                    fprintf(fileID2,'\n');
                end
            end
        end
        
    end
    fclose(fileID1);
    fclose(fileID2);
end

function [lambda_new, pitch_new, Cp_new, Ct_new] = upscale(lambda, pitch, Cp, Ct, n1, n2)
    lambda_new = linspace(min(lambda),max(lambda),n1);
    pitch_new = linspace(min(pitch),max(pitch),n2)';
    Cp_new = interp2(pitch, lambda, Cp', pitch_new, lambda_new)';
    Ct_new = interp2(pitch, lambda, Ct', pitch_new, lambda_new)';
end