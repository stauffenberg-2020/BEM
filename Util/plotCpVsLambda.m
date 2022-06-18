function [CM,cc] = plotCpVsLambda(pitch, lambda, Cp)

    [x,y]=find(ismember(Cp,max(Cp(:))));
    idx_to_plot = unique([flip(x:-8:1) x:6:length(pitch)]); % Making sure max Cp is plotted
    j=1;
    for i=[idx_to_plot]
        [max_Cp(j),idx(j)] = max(Cp(i,:));
        j=j+1;
    end
    [CM,cc] = plotCpLambdaContour(lambda, pitch(idx_to_plot)', (Cp(idx_to_plot,:))');
    hold on
    ids_to_plot = reduce(max_Cp);
    plot(lambda(idx(ids_to_plot)),max_Cp(ids_to_plot),'.-r','HandleVisibility','off');
    plot(lambda(y),max(Cp(:)),'r-p')
    txt = sprintf('  Cp_m_a_x = %0.3f',max(Cp(:)));
    text(lambda(y),max(Cp(:)),txt,'Color','red');
end

%% Supporting function
function [CM,cc] = plotCpLambdaContour(inp_x,inp_z,inp_y)
    
    % inp_x: 1xM-vector with x-coordinates
    % inp_z: 1xN-vector with z-coordinates (the 'text-Info')
    % inp_y: MxN-matrix with datapoints

    [X,Y]=meshgrid(inp_x,inp_z);
    [X_inp,Y_inp]=meshgrid(inp_x,inp_z);
    Z=interp2(X_inp,Y_inp,inp_y',X,Y);
    [CM, cc] = contour(X,Z,Y,inp_z,'ShowText','on');
    ylim([0 0.593]);
    grid on
    xlabel('Lambda')
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