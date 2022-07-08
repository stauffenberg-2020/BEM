function [C,h] = plotPitchVsLambda(pitch, lambda, Cp)
    [C,h] = contourf(lambda',pitch,Cp,500,'LineStyle','--');
    hold on
    xlabel('Lambda (TSR)');
    ylabel('Pitch (deg)');

    Lvls = h.LevelList;
    h.LevelList = Lvls(Lvls >= 0); 
    h.LevelList = round(h.LevelList,3);
    clabel(C,h);

    [x,y]=find(ismember(Cp,max(Cp(:))));
    plot(lambda(y),pitch(x),'r-p')
    txt = sprintf('  Cp_m_a_x = %0.3f',max(Cp(:)));
    text(lambda(y),pitch(x),txt,'Color','red');
    line([lambda(1),lambda(y)],[pitch(x),pitch(x)],'Color','red','LineStyle','--')
    line([lambda(y),lambda(y)],[pitch(1),pitch(x)],'Color','red','LineStyle','--')
    hcb = colorbar;
    hcb.Title.String = 'Cp';

    for i=1:length(lambda)
        [~,I(i)] = max(Cp(:,i));
    end
    plot(lambda,pitch(I),'.-r','LineWidth',0.5);
    
    
end


