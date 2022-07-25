function h = plot_over_disc(BEM, op_pt_id ,param)
    h = polar([0 2*pi], [0 BEM.r(end)]);
    view(90,-270);
    delete(h); % check for alternates
    
%     polarplot([0 2*pi], [0 BEM.r(end)]);
%     ax = gca;
%     ax.ThetaZeroLocation = 'bottom';

    hold on;

    bld_no = 1; % Plotting parameters extracted from blade 1
    b = {'b1','b2','b3'};
    
    fn = fieldnames(BEM.(b{bld_no}));
    azim_flag = startsWith(fn,'az_');
    azimuths = fn(azim_flag);
    azims = str2double(extractAfter(azimuths,'az_'));
    
    r = BEM.r;
    
    for i=1:length(azims)
        az = ['az_',num2str(azims(i),'%03d')];
        param_azi(:,i) = BEM.(b{bld_no}).(az).(param)(:,op_pt_id);
    end
   
    [TH,R] = meshgrid(azims,r);
    [X,Y] = pol2cart(deg2rad(TH),R);
    surf(X,Y,param_azi,'EdgeColor','flat','FaceColor','interp');
    hcb = colorbar;
    title(hcb,strrep(param,'_','-'));
end