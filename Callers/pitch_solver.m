% Thrust Limiter Investigations
clearvars;clc;
op_data = readmatrix('G:\BEM\BEM\Data\NREL5MWRefTurb_v50\data\operational_data.opt','FileType','text');
ae = readmatrix('G:\BEM\BEM\Data\NREL5MWRefTurb_v50\data\NREL_5MW_ae.txt');
pc = 'G:\BEM\BEM\Data\NREL5MWRefTurb_v50\data\NREL_5MW_pc.txt';

N=3;
rho = 1.225; % Air Density [Kg/m^3]
op_pts.wsp = op_data(:,1);
op_pts.pitch = op_data(:,2); % Operational set points, Collective Pitch [deg]
op_pts.rpm = op_data(:,3); % Operational set points, LSS Rotor speed [rpm]

sec_r = ae(2:size(ae,1),1);
sec_C = ae(2:size(ae,1),3);
sec_t_C = ae(2:size(ae,1),5);
sec_twist = [-13.3,-13.3,-13.3,-13.3,-11.48,-10.16,-9.011,-7.795,-6.544,-5.361,-4.188,-3.125,-2.319,-1.526,-0.863,-0.37,-0.106, 0]'; % to be read from *.htc file


BLD = organize_blade_data(sec_r, sec_C, sec_t_C, sec_twist, pc);


AeroFlags.induction = 0; % 0 or 1, 0 = Induction Off, 1 = Induction On
AeroFlags.tip_loss = 0; % 0 or 1, 0 = Prandtl's tip loss correction Off, 1 = Prandtl's tip loss correction On
AeroFlags.highCT = 0; % 0 or 1 or 2, 0 = Off, 1 = As per HANSEN eqn. 6.38, 2 = As per HANSEN eqn. 6.37

%% Core BEM calculations

[output_details,output,BEM] = core_bem(rho, N, op_pts, BLD, AeroFlags);

%% Finding solution for target Root Flap Moment
target = max(output(:,3))*[0.96 0.97 0.98 0.99]; % Maximum root flap load target [kNm]

for i=1:length(target)
    tuning_ids{:,i}=find(output(:,3)>target(i)); % finding indices of operating points that exceed the target

    if ~isempty(tuning_ids{i})
        % Better method by using optimization toolbox
%         fun=@(x)solve_for_pitch(rho, op_pts.wsp(tuning_ids{i}), x, op_pts.rpm(tuning_ids{i}), BLD, AeroFlags)-target(i);
%         tuned_pitch = lsqnonlin(fun,op_pts.pitch(tuning_ids{i})); % Solving the function with a initial value

        % Basic method within base MATLAB
        fun=@(x)sum((solve_for_pitch(rho, op_pts.wsp(tuning_ids{i}), x, op_pts.rpm(tuning_ids{i}), BLD, AeroFlags)-target(i)).^2);
        tuned_pitch = fminsearch (fun,op_pts.pitch(tuning_ids{i})); 

        % Updating the output
        op_pts.pitch_tuned{i} = op_pts.pitch; % Creating a copy of the original output (for further updation)
        op_pts.pitch_tuned{i}(tuning_ids{i}) = tuned_pitch; % Tuned values updated
        
        op_pts_tuned.wsp = op_pts.wsp(tuning_ids{i});
        op_pts_tuned.pitch = tuned_pitch;
        op_pts_tuned.rpm = op_pts.rpm(tuning_ids{i});

        [~,output_tune] = core_bem(rho, N,op_pts_tuned, BLD, AeroFlags);
        output_tuned{i} = output; % Creating a copy of the original output (for further updation)
        output_tuned{i}(tuning_ids{i},:) = output_tune; % Updated output
    end
end
%% Plotting

figure('Renderer', 'painters', 'Position', [10 10 1000 1000])
subplot(2,2,1)
plot(op_pts.wsp,op_pts.pitch,'k*-','DisplayName','TL Off')
grid on
hold on
for i=1:length(target)
    plot(op_pts.wsp,op_pts.pitch_tuned{i},'*-','color',Plot_cols(i),'DisplayName',num2str(round(target(i))));
end
xlabel('Wind Speed (m/s)');
ylabel('Pitch (deg)');
leg = legend('location','best');
title(leg,'-Mx11h target')

subplot(2,2,2)
plot(op_pts.wsp,output(:,3),'k*-','DisplayName','TL Off')
grid on
hold on
for i=1:length(target)
    plot(op_pts.wsp,output_tuned{i}(:,3),'*-','color',Plot_cols(i),'DisplayName',num2str(round(target(i))));
end
xlabel('Wind Speed (m/s)');
ylabel('-Mx11h (kNm)');
leg = legend('location','best');
title(leg,'-Mx11h target')


P_rot_TLOff = output(:,5).*(2*pi.*op_pts.rpm/60); % Rotor power calculated from Maero
M=8;
k=2.5;
[AEP_MWh_TLOff, ~] = AEP(op_pts.wsp,P_rot_TLOff,M,k);

hfig=subplot(2,2,3);
yyaxis left
plot(op_pts.wsp,P_rot_TLOff,'k*-','DisplayName','TL Off');
hold on
grid on
for i=1:length(target)
    P_rot_tuned(:,i) = output_tuned{i}(:,5).*(2*pi.*op_pts.rpm/60); % Rotor power calculated from Maero
    [AEP_MWh_tuned(i), ~] = AEP(op_pts.wsp,P_rot_tuned(:,i),M,k);
    plot(op_pts.wsp,P_rot_tuned(:,i),'*-','color',Plot_cols(i),'DisplayName',num2str(round(target(i))));
end
xlabel('Wind Speed (m/s)');
ylabel('Rotor Power (kW)');
leg = legend('location','se');
title(leg,'-Mx11h target')
title(sprintf('M = %0.1f, k = %0.2f',M,k))

yyaxis right
weib_pdf = wblpdf(op_pts.wsp,M/(gamma(1+1/k)),k);
plot(op_pts.wsp, weib_pdf,'k--','handlevisibility','off');
hfig.YAxis(1).Color = 'k';
hfig.YAxis(2).Visible = 'off';


% AEP drop calculation
AEP_drop = ((AEP_MWh_TLOff-AEP_MWh_tuned(:))/AEP_MWh_TLOff)*100;
load_drop = ((max(output(:,3))-target)./max(output(:,3)))*100;
subplot(2,2,4)
plot(AEP_drop,load_drop,'k--')
grid on
hold on
for i=1:length(AEP_drop)
    plot(AEP_drop(i),load_drop(i),'*','color',Plot_cols(i))
end
xlabel('AEP Drop in %');
ylabel('Flap Load drop in %');
P = polyfit(AEP_drop,load_drop,1);
txt = {sprintf('Slope = %0.2f',P(1));sprintf('Intercept = %0.2f',P(2))};
xl=xlim;
yl=ylim;
text(mean(xl)+0.01/mean(xl),mean(ylim),txt)
title(sprintf('M = %0.1f, k = %0.2f',M,k))

%% Supporting functions

function out = solve_for_pitch(rho, wsp, x, rpm, BLD, AeroFlags)
    N=3;
    op_pts.wsp = wsp;
    op_pts.pitch = x; % Operational set points, Collective Pitch [deg]
    op_pts.rpm = rpm; % Operational set points, LSS Rotor speed [rpm]
    [~, results] = core_bem(rho, N, op_pts, BLD, AeroFlags);
    out = results(:,3); % Choosing the flap load as target output
end




function [AEP_MWh, Energy] = AEP(WSP,P,M,k)
    A                  = M/(gamma(1+1/k));
    CDF                = 1-exp(-(WSP/A).^k);
    N_CDF              = length(CDF);
    for ii = 1:N_CDF-1
        Time(ii)        = CDF(ii+1)-CDF(ii);
        MeanP(ii)       = (P(ii+1)+P(ii))/2;
    end
    Energy             = Time.*MeanP;
    AEP_MWh            = sum(Energy)*8760/1000;
end


function out = Plot_cols(i)
    col = [0    0.4470    0.7410;  % Blue
        0.8500    0.3250    0.0980; % Red
        0.9290    0.6940    0.1250; % Yellow
        0.4940    0.1840    0.5560; % Purple
        0.4660    0.6740    0.1880; % Green
        0.3010    0.7450    0.9330; % Light Blue
        0.6350    0.0780    0.1840];% Dark Red
    
    i = mod(i-1,length(col(:,1)))+1;
    out = col(i,:);
end







