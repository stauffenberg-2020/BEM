function [pitch, region] = core_pitch_ctr(wsp, rpm, rho, R, pitch_range, lambda_range, Cp_matrix, min_rpm, rtd_rpm, Prat, eff)

    % Simple collective pitch controller based on Cp data
    % 
    % Outputs
    % Pitch [deg] is the collective pitch for each supplied wsp
    % region number of the controller (1 to 4) for each supplied wind speed
    
    I = zeros(length(lambda_range),1);
    for i=1:length(lambda_range)
        [~,I(i)] = max(Cp_matrix(:,i)); % Index of pitch for maximum Cp for each Lambda
    end
    OTC_X = lambda_range;
    OTC_Y = pitch_range(I); % Pitch corresponding to max Cp for for each lambda, OTC_Y
    
    lambda_wsp = (rpm'.*2*pi/60).*R./wsp'; % (local) Lambda for the given wsp wind speeds
    region = zeros(length(rpm),1);
    pitch = zeros(length(rpm),1);
    CP = zeros(length(rpm),1);
    
    for i=1:length(rpm)
        % Calculations for region 3 check
        rtd_rot_pow = el_pow_to_rot_pow(Prat, eff);
        % Calculating Cp value corresponding to lambda(i) & pitch(i). pitch(i) is from regular OTC table to maximize Cp
        Cp = interp2(lambda_range, pitch_range, Cp_matrix, lambda_wsp(i), interp1(OTC_X, OTC_Y,lambda_wsp(i),'spline')); %  Check this
        rot_pow = (0.5*rho*pi*R*R.*wsp(i).^3).*Cp*(1/1000);
        
        if rpm(i) == rtd_rpm && rot_pow < rtd_rot_pow
            region(i,1) = 3; % Constant rpm region, but not yet constant power
        elseif rpm(i) == rtd_rpm
            region(i,1) = 4; % Constant rpm & constant power region
        elseif rpm(i) == min_rpm
            region(i,1) = 1; % Constant rpm region
        else 
            region(i,1) = 2; % Varying rpm region
        end
        
        
        switch region(i) % The 4 different regions of controller
            case 1 % Constant rpm region
                pitch(i,1) = interp1(OTC_X, OTC_Y, lambda_wsp(i),'spline'); % Check this
            case 2 % Varying rpm region
                pitch(i,1) = interp1(OTC_X, OTC_Y, lambda_wsp(i));
            case 3 % Constant rpm region, but not yet constant power, So maintaining pitch corresponding to max Cp
                pitch(i,1) = pitch(i-1,1); % Check this
%                 pitch(i,1) = interp1(OTC_X,OTC_Y,lambda_wsp(i)); % Check this
            case 4 % Constant rpm & constant power region
                CP(i) = el_pow_to_rot_pow(Prat, eff)*1000./(0.5*rho*pi*R*R.*wsp(i).^3); % Full load pitch calculation
                CPvsPit = interp2(lambda_range,pitch_range,Cp_matrix,repmat(lambda_wsp(i),length(pitch_range),1),pitch_range); % CP as a function of pitch (CP vs pitch) for a given lambda
                pit = intersect(CPvsPit, pitch_range, CP(i));
                pitch(i,1) = max(pit); % Check max on min or end
        end
     end
    
end

%% Supporting function(s)
function P_rot = el_pow_to_rot_pow(Pel, eff)
    P_rot = Pel/eff;
end






