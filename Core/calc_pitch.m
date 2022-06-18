function [pitch, region] = calc_pitch(wsp, rpm, rho, R, pitch, lambda, Cp_matrix, min_rpm, rtd_rpm, Prat)
    % Simple collective pitch controller based on Cp data
    % 
    % Outputs
    % Pitch [deg] is the collective pitch for each supplied wsp
    % region number of the controller (1 to 4) for each supplied wind speed
    
    for i=1:length(lambda)
        [~,I(i)] = max(Cp_matrix(:,i)); % Index of pitch for maximum Cp for each Lambda
    end
    OTC_X = lambda;
    OTC_Y = pitch(I); % Pitch corresponding to max Cp for for each lambda, OTC_Y
    
    Cp_data.pitch = pitch;
    Cp_data.lambda = lambda;
    Cp_data.Cp= Cp_matrix;
    
    lambda_wsp = (rpm'.*2*pi/60).*R./wsp'; % (local) Lambda for the given wsp wind speeds
    region = zeros(length(rpm),1);
    pitch = zeros(length(rpm),1);
    
    for i=1:length(rpm)
        % Calculations for region 3 check
        rtd_rot_pow = el_pow_to_rot_pow(Prat);
        % Calculating Cp value corresponding to lambda(i) & pitch(i). pitch(i) is from regular OTC table to maximize Cp
        Cp = interp2(Cp_data.lambda, Cp_data.pitch, Cp_data.Cp, lambda_wsp(i), interp1(OTC_X, OTC_Y,lambda_wsp(i),'spline')); %  Check this
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
%                 pitch(i,1) = interp1(OTC_X,OTC_Y,lambda(i)); % Check this
            case 4 % Constant rpm & constant power region
                CP(i) = el_pow_to_rot_pow(Prat)*1000./(0.5*rho*pi*R*R.*wsp(i).^3); % Full load pitch calculation
                CPvsPit = interp2(Cp_data.lambda,Cp_data.pitch,Cp_data.Cp,lambda_wsp(i),Cp_data.pitch); % CP as a function of pitch (CP vs pitch) for a given lambda
                [~,pit] = intersections(CPvsPit, Cp_data.pitch,[CP(i) CP(i)],[min(Cp_data.pitch) max(Cp_data.pitch)]);
                pitch(i,1) = min(pit);
        end
     end
    
end

%% Supporting function(s)
function P_rot = el_pow_to_rot_pow(Pel)
    P_rot = Pel/0.944;
end






