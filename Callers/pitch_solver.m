% Flap Load Limiting Investigations
clearvars;clc;

filename = 'NREL_5MW.txt';
[General, op_pts, BLD, ~] = read_turbine_file(filename);
BLD.preflap = zeros(length(BLD.r),1);
op_pts.pitch(9:end) = [3.83,6.60,8.7,10.45,12.06,13.54,14.92,16.23,17.47,18.70,19.94,21.18,22.35,23.47]'; % Values from https://www.nrel.gov/docs/fy09osti/38060.pdf

[output_details, output, BEM] = core_bem(General, op_pts, BLD);

%% Finding solution for target Root Flap Moment
target = max(output(:,3))*[0.8 0.85 0.9 0.95]; % Maximum root flap load target [kNm]

for i=1:length(target)
    tuning_ids{:,i}=find(output(:,3)>target(i)); % finding indices of operating points that exceed the target

    if ~isempty(tuning_ids{i})
        
        tuning_pts.wsp = op_pts.wsp(tuning_ids{i});
        tuning_pts.pitch = op_pts.pitch(tuning_ids{i});
        tuning_pts.rpm = op_pts.rpm(tuning_ids{i});

        
        % Better method by using optimization toolbox
        fun=@(x)solve_for_pitch(x, General, tuning_pts, BLD)-target(i);
        tuned_pitch = lsqnonlin(fun,tuning_pts.pitch); % Solving the function with a initial value

        % Basic method within base MATLAB
%         fun=@(x)sum((solve_for_pitch(x, General, tuning_pts, BLD)-target(i)).^2);
%         tuned_pitch = fminsearch (fun,tuning_pts.pitch); 

        % Updating the output
        op_pts.pitch_tuned{i} = op_pts.pitch; % Creating a copy of the original output (for further updation) 
        op_pts.pitch_tuned{i}(tuning_ids{i}) = tuned_pitch; % Tuned values updated         
        
        op_pts_tuned.wsp = op_pts.wsp(tuning_ids{i});
        op_pts_tuned.pitch = tuned_pitch;
        op_pts_tuned.rpm = op_pts.rpm(tuning_ids{i});
        
        [~,output_tune] = core_bem(General, op_pts_tuned, BLD);
        output_tuned{i} = output; % Creating a copy of the original output (for further updation)
        output_tuned{i}(tuning_ids{i},:) = output_tune; % Updated output
    end
    clear tuning_pts
end
%% Plotting

MP = get(0,'MonitorPositions');
if size(MP,1) == 1
    SS = MP(1,:);
else
    SS = MP(2,:); % To have the results displayed better in bigger screen
end
figure('position',[SS(1)+SS(3)*0.15, SS(2)+SS(4)*0.15, SS(3)*0.7, SS(4)*0.7]);
subplot(2,2,1)
plot(op_pts.wsp,op_pts.pitch,'k*-','DisplayName','Base')
grid on
hold on
for i=1:length(target)
    plot(op_pts.wsp,op_pts.pitch_tuned{i},'*-','color',plot_cols(i),'DisplayName',[num2str(round(target(i)/max(output(:,3)),2)) ' %']);
end
xlabel('Wind Speed (m/s)');
ylabel('Pitch (deg)');
leg = legend('location','nw');
title(leg,{'Max Flap' ,'Load target'})

subplot(2,2,2)
plot(op_pts.wsp,output(:,3),'k*-','DisplayName','Base')
grid on
hold on
for i=1:length(target)
    plot(op_pts.wsp,output_tuned{i}(:,3),'*-','color',plot_cols(i),'DisplayName',[num2str(round(target(i)/max(output(:,3)),2)) ' %']);
end
xlabel('Wind Speed (m/s)');
ylabel('Root flap load (kNm)');
% leg = legend('location','best');
% title(leg,'Load target')


P_rot_TLOff = output(:,5).*(2*pi.*op_pts.rpm/60); % Rotor power calculated from Maero
M=8;
k=2.5;
[AEP_MWh_TLOff, ~] = AEP(op_pts.wsp,P_rot_TLOff,M,k);

hfig=subplot(2,2,3);
yyaxis left
plot(op_pts.wsp,P_rot_TLOff,'k*-','DisplayName','Base');
hold on
grid on
for i=1:length(target)
    P_rot_tuned(:,i) = output_tuned{i}(:,5).*(2*pi.*op_pts.rpm/60); % Rotor power calculated from Maero
    [AEP_MWh_tuned(i), ~] = AEP(op_pts.wsp,P_rot_tuned(:,i),M,k);
    plot(op_pts.wsp,P_rot_tuned(:,i),'*-','color',plot_cols(i),'DisplayName',[num2str(round(target(i)/max(output(:,3)),2)) ' %']);
end
xlabel('Wind Speed (m/s)');
ylabel('Rotor Power (kW)');
% leg = legend('location','se');
% title(leg,'Load target')
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
    plot(AEP_drop(i),load_drop(i),'*','color',plot_cols(i))
end
xlabel('AEP Drop in % (w.r.t. Base)');
ylabel('Flap Load drop in %');
P = polyfit(AEP_drop,load_drop,1);
txt = {sprintf('Slope = %0.2f',P(1));sprintf('Intercept = %0.2f',P(2))};
xl=xlim;
yl=ylim;
% text(mean(xl),mean(ylim),txt)
title(sprintf('M = %0.1f, k = %0.2f',M,k))

%% Supporting functions

function flap_loads = solve_for_pitch(x, General, op_pts, BLD)
    op_pts.pitch = x; % Operational set points, Collective Pitch [deg]
    [~, results,~] = core_bem(General, op_pts, BLD);
    flap_loads = results(:,3); % Choosing the flap load as target output
end















