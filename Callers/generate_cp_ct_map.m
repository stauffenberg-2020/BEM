% Script for generating Cp (and Ct) map

% Caller script
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

lambda = 2:1:20; % Lambda (TSR) range
op_pts.pitch = (-5:0.5:25)'; % Pitch range
wsp = 8; % Just a temp variable for deriving rpm corresponding to lambda

Cp = zeros(length(op_pts.pitch),length(lambda));
Ct = zeros(length(op_pts.pitch),length(lambda));
for i=1:length(op_pts.pitch)
    map_pts.wsp(1:length(lambda),1) = wsp;
    map_pts.pitch(1:length(lambda),1) = op_pts.pitch(i);
    map_pts.rpm = (lambda'*wsp/BLD.r(end)).*(30/pi);
    [~, output, ~] = core_bem(General, map_pts, BLD);
    Cp(i,:) = output(:,7);
    Ct(i,:) = output(:,9);
end
%%
path = '.\Data\';
write_cp_ct(lambda,op_pts.pitch,Cp,Ct,path);
%% Plotting

close all;
MP = get(0,'MonitorPositions');
if size(MP,1) == 1
    SS = MP(1,:);
else
    SS = MP(2,:); % To have the results better displayed in bigger screen
end
figure('position',[SS(1)+SS(3)*0.05, SS(2)+SS(4)*0.2, SS(3)*0.9, SS(4)*0.6])

subplot(1,2,1)
[C,h] = contourf(lambda',op_pts.pitch,Cp,500,'LineStyle','--');
hold on
xlabel('Lambda (TSR)');
ylabel('Pitch (deg)');

Lvls = h.LevelList;
h.LevelList = Lvls(Lvls >= 0); 
h.LevelList = round(h.LevelList,3);
clabel(C,h);

[x,y]=find(ismember(Cp,max(Cp(:))));
plot(lambda(y),op_pts.pitch(x),'r-p')
txt = sprintf('  Cp_m_a_x = %0.3f',max(Cp(:)));
text(lambda(y),op_pts.pitch(x),txt,'Color','red');
line([lambda(1),lambda(y)],[op_pts.pitch(x),op_pts.pitch(x)],'Color','red','LineStyle','--')
line([lambda(y),lambda(y)],[op_pts.pitch(1),op_pts.pitch(x)],'Color','red','LineStyle','--')
hcb = colorbar;
hcb.Title.String = 'Cp';

for i=1:length(lambda)
    [M(i),I(i)] = max(Cp(:,i));
end
plot(lambda,op_pts.pitch(I),'.-r','LineWidth',0.5);

subplot(1,2,2)
idx_to_plot = unique([flip(x:-4:1) x:4:length(op_pts.pitch)]); % Making sure max Cp is plotted
j=1;
for i=[idx_to_plot]
    [max_Cp(j),idx(j)] = max(Cp(i,:));
    j=j+1;
end
[CM,cc] = plotCpVsLambda(lambda,op_pts.pitch(idx_to_plot)',(Cp(idx_to_plot,:))');
hold on
ids_to_plot = reduce(max_Cp);
plot(lambda(idx(ids_to_plot)),max_Cp(ids_to_plot),'.-r','HandleVisibility','off');
plot(lambda(y),max(Cp(:)),'r-p')
text(lambda(y),max(Cp(:)),txt,'Color','red');

%% Supporting functions 
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

function plot_ids = reduce(Cp)
    Cp_red = Cp;
    X = Cp(length(Cp));
    for k = length(Cp)-1:-1:1
      if Cp(k) < X
        Cp_red(k) = X;  % Or: NaN
      else
        X = Cp(k);
      end
    end
    plot_ids = (Cp_red == Cp);
end

function write_cp_ct(lambda,pitch,Cp,Ct,path)
    %% Writing the outputs into a text file
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
                    fprintf(fileID1,'%8.3f \t',lambda(j-1));
                    fprintf(fileID2,'%8.3f \t',lambda(j-1));
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
                    fprintf(fileID1,'%8.3f \t',Cp(i-1,j-1));
                    fprintf(fileID2,'%8.3f \t',Ct(i-1,j-1));
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