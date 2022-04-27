function [output_details, output, BEM] = core_bem(rho, N, wsp, pitch, rpm, sec_r, sec_C, sec_t_C, sec_twist, profile_t_C, profile_AoA, profile_cL, profile_cD , induction_flag, tip_loss_flag, glauert_corr_flag)
    % Core steady state BEM calculations
    % Inputs
    % 
    % rho [Kg/m^3] is the air density
    % N is the number of blades
    % wsp [m/s] is the wind speed (can be a n x 1 vector with multiple entries)
    % Pitch [deg] is the collective pitch (can be a n x 1 vector with multiple entries, But should have 1:1 correspondance with wsp)
    % rpm [rpm] is the (LSS) rotor speed in rpm (can be a nx1 vector with multiple entries, But should have 1:1 correspondance with wsp)
    % sec_r [m] is is the blade section table (m x 1 vector)
    % sec_C [m] is the Chord corresponding to each of the sec_r (m x 1 vector)
    % sec_t_C is the thickness to chord ratio corresponding to each of the sec_r (m x 1 vector)
    % profile_t_C is the K x 1 thickneses available
    % profile_AoA is the L x 1 angles of attack
    % profile_cL is L x K cL values corresponding to each angles of attack and t_C
    % profile_cD is L x K cD values corresponding to each angles of attack and t_C 
    % induction_flag for both axial and tangential induction, 0 = Induction off, 1 = Induction on
    % tip_loss_flag is Prandtl's Tip correction flag, 0 = correction off, 1 = correction on
    % glauert_corr_flag is Galuert's High CT correction, 0 = Off, 1 = As per HANSEN eqn. 6.38, 2 = As per HANSEN eqn. 6.37
    
    % Input data organising, Making sure values are in column vector
    wsp = wsp(:); 
    pitch = pitch(:);
    rpm = rpm(:);
    sec_r = sec_r(:);
    sec_C = sec_C(:);
    sec_t_C = sec_t_C(:);
    sec_twist = sec_twist(:);
    profile_t_C = profile_t_C(:);
    profile_AoA = profile_AoA(:);

    
    BEM.r = sec_r; % Blade stations
    BEM.C = sec_C; % Chord distribution
    BEM.t_C = sec_t_C; % t/C distribution
    BEM.AeroTwist = sec_twist; % Aerodynamic (mold) Twist distribution [deg]
    
    BEM.Omega_r = (rpm'*2*pi/60).*BEM.r; % Velocity seen by each of the section due to rotation [m/s]
    BEM.aA = zeros(length(BEM.r),size(wsp,1)); % Initializing axial induction factors 
    BEM.aT = zeros(length(BEM.r),size(wsp,1)); % Initializing tangential induction factors 
    BEM.V_Res = sqrt(wsp'.^2+BEM.Omega_r.^2); % Resultant velocity seen by each of the section [m/s]
    BEM.phi = rad2deg(atan2((wsp'.*(1-BEM.aA)),(BEM.Omega_r.*(1+BEM.aT)))); % Calculating flow angle PHI [deg]
    BEM.alpha = BEM.phi-BEM.AeroTwist-pitch'; % Angle of attack seen by each section [deg]
    [BEM.CL, BEM.CD] = CL_CD_vs_alpha(BEM.t_C, profile_t_C, profile_AoA, profile_cL, profile_cD, BEM.alpha); % CL, CD coefficients
    BEM.sig = BEM.C*N./(2*pi.*BEM.r); % Solidity ratio
    BEM.CN = BEM.CL.*cos(deg2rad(BEM.phi))+BEM.CD.*sin(deg2rad(BEM.phi)); % Sectional Normal force coefficient
    BEM.CT = BEM.CL.*sin(deg2rad(BEM.phi))-BEM.CD.*cos(deg2rad(BEM.phi)); % Sectional Tangential force coefficient
    try
        if induction_flag
            for i=1:length(BEM.r) % Loop for different sections of the blade
                for j=1:size(wsp,1) % Loop for different operating set points
                    BEM.aA_new = 1./((4.*sin(deg2rad(BEM.phi)).*sin(deg2rad(BEM.phi))./(BEM.sig.*BEM.CN))+1); % Calculating new axial induction factors
                    BEM.aT_new = 1./((4.*sin(deg2rad(BEM.phi)).*cos(deg2rad(BEM.phi))./(BEM.sig.*BEM.CT))-1); % Calculating new axial induction factors
                    iter=1;
                    while (abs(BEM.aA_new(i,j) - BEM.aA(i,j)) > 0.0001 || abs(BEM.aT_new(i,j) - BEM.aT(i,j)) > 0.0001) && iter < 1000 % Induction iteration
                        if BEM.aA_new(i,j) > 1.5, BEM.aA(i,j) = 1.5; elseif BEM.aA_new(i,j) < -1, BEM.aA(i,j) = -1; else, BEM.aA(i,j) = BEM.aA_new(i,j);end
                        if BEM.aT_new(i,j) > 1,   BEM.aT(i,j) = 1;   elseif BEM.aT_new(i,j) < -1, BEM.aT(i,j) = -1; else, BEM.aT(i,j) = BEM.aT_new(i,j);end
                        BEM.phi(i,j) = rad2deg(atan2((wsp(j)'.*(1-BEM.aA(i,j))),(BEM.Omega_r(i,j).*(1+BEM.aT(i,j))))); % Calculating flow angle PHI [deg]
                        if tip_loss_flag
                            f(i,j) = (N/2).*(BEM.r(end)-BEM.r(i))./(BEM.r(i).*sin(deg2rad(BEM.phi(i,j))));
                            BEM.F(i,j)=(2/pi).*acos(exp(-f(i,j)));
                        else
                            BEM.F(i,j)=1;
                        end
                        BEM.alpha(i,j) = BEM.phi(i,j)-BEM.AeroTwist(i)-pitch(j)'; % Angle of attack seen by each section [deg]
                        if BEM.alpha(i,j) < -180 || BEM.alpha(i,j) > 180
                            BEM.alpha(i,j) = mod(BEM.alpha(i,j), 180); % Check here once
                        end
                        [BEM.CL(i,j), BEM.CD(i,j)] = CL_CD_vs_alpha(BEM.t_C(i), profile_t_C, profile_AoA, profile_cL, profile_cD, BEM.alpha(i,j)); % CL, CD coefficients
                        BEM.CN(i,j) = BEM.CL(i,j).*cos(deg2rad(BEM.phi(i,j)))+BEM.CD(i,j).*sin(deg2rad(BEM.phi(i,j))); % Sectional Normal force coefficient (Thrust direction force coefficient)
                        BEM.CT(i,j) = BEM.CL(i,j).*sin(deg2rad(BEM.phi(i,j)))-BEM.CD(i,j).*cos(deg2rad(BEM.phi(i,j))); % Sectional Tangential force coefficient
                        if glauert_corr_flag
                            method = glauert_corr_flag; % Method = 1 or 2. Method 2 is similar to Pyro. Method 1 is similar to VTS
                            BEM.aA_new(i,j) = calc_ind_using_high_CT_approx(method, BEM.aA(i,j), BEM.F(i,j), BEM.phi(i,j), BEM.sig(i), BEM.CN(i,j));
                        else
                            BEM.aA_new(i,j) = 1./((4.*BEM.F(i,j).*sin(deg2rad(BEM.phi(i,j))).*sin(deg2rad(BEM.phi(i,j)))./(BEM.sig(i).*BEM.CN(i,j)))+1); % Calculating new axial induction factors
                        end
                        BEM.aT_new(i,j) = 1./((4.*BEM.F(i,j).*sin(deg2rad(BEM.phi(i,j))).*cos(deg2rad(BEM.phi(i,j)))./(BEM.sig(i).*BEM.CT(i,j)))-1); % Calculating new axial induction factors
                        iter=iter+1;
                    end
                    BEM.iter_count(i,j)=iter; % Induction iteration counter
                end
            end
        end
        disp('BEM calculations successfully completed')
    catch ME
        fprintf('ERROR in blade station %0.3f(m) at wind speed = %0.2f(m/s) during induction iteration %d \n', round(sec_r(i),2), wsp(j), iter)
        if BEM.phi(i,j) < 0
            fprintf('Possible details : Axial induction of %0.2f & Tangential Induction = %0.2f results in negative flow angle of %0.2f \n', BEM.aA(i,j), BEM.aT(i,j), BEM.phi(i,j))
        end
        fprintf('Debugging details : \n')
        for i = 1:numel(ME.stack)
            disp(ME.stack(i))
        end
        fprintf('Warning : Outputs calculated only till error occurrence \n')
    end
    
    %% Sectional Calculations
    A = pi.*BEM.r(end).^2; % Actual Rotor area
    Q_L = 0.5*rho.*BEM.V_Res.*BEM.V_Res; % Local Dynamic pressure for each station
    Q_G = 0.5.*rho.*wsp'.*wsp'; % Dynamic pressure
    
    BEM.FlapdF = Q_L.*BEM.C.*BEM.CN; % Local Thrust/Normal force [N], Check this
    BEM.FlapdM = BEM.FlapdF.*(BEM.r-BEM.r(1)); % Local Flap Moment at Hub end or Blade root [Nm]
    BEM.EdgedF = Q_L.*BEM.C.*BEM.CT; % Local Edge force [N]
    BEM.EdgedM = BEM.EdgedF.*(BEM.r-BEM.r(1)); % Local Edge Moment at Hub end or Blade root [Nm]
    BEM.EdgedMC = BEM.EdgedF.*BEM.r; % Local Edge Moment at rotor center [Nm]
    
    BEM.dT = N.*BEM.FlapdF;
    BEM.Ct = BEM.dT./(Q_G.*2.*pi.*BEM.r); % Local Thrust Coefficient from Local Thrust
    

    
    %% Output calculations
    
    output(:,1) = (1/1000)*trapz(BEM.r(1:length(BEM.r)-1,:),BEM.FlapdF(1:length(BEM.r)-1,:)); % Flap Force [KN], % Ignoring the influence of last section similar to VTS
    output(:,2) = (1/1000)*trapz(BEM.r(1:length(BEM.r)-1,:),BEM.EdgedF(1:length(BEM.r)-1,:)); % Edge Force [KN]
    output(:,3) = (1/1000)*trapz(BEM.r(1:length(BEM.r)-1,:),BEM.FlapdM(1:length(BEM.r)-1,:)); % Flap Moment [KNm]
    output(:,4) = (1/1000)*trapz(BEM.r(1:length(BEM.r)-1,:),BEM.EdgedM(1:length(BEM.r)-1,:)); % Edge Moment [KNm]
    output(:,5) = (1/1000)*trapz(BEM.r(1:length(BEM.r)-1,:),BEM.EdgedMC(1:length(BEM.r)-1,:))*N; % Aerodynamic Torque [KNm]
    output(:,6) = output(:,5).*(2*pi.*rpm/60); % Rotor Power [KW]
    output(:,7) = (output(:,6).*1000)./(Q_G'.*A.*wsp); % Rotor averaged Cp calculated from rotor power
    output(:,8) = output(:,1)*N; % Aerodynamic Thrust [KN]
    output(:,9) = (output(:,8).*1000)./(Q_G'.*A); % Rotor averaged Ct Calculated from rotor thrust
    
    output_details = {'Fy11h', 'Fx11h', '-Mx11h', 'My11h', 'Maero', 'P_rot', 'Cp', 'Fthr', 'Ct';
                            'Flap Force', 'Edge Force', 'Flap Moment', 'Edge Moment', 'Aerodynamic Torque', 'Rotor Aerodynamic Power', 'Rotor Cp','Rotor Aerodynamic Thrust','Rotor Ct';
                            '[kN]','[kN]','[kNm]','[kNm]','[kNm]','[kW]','[-]','[kN]','[-]'};
    


end

%% Supporting functions

function index = first_higher(A, target)
    b = (find(A > target));
    index = b(1);
end

function index = first_lower(A, target)
    b = (find(A < target));
    index = b(end);
end

function [CL, CD] = CL_CD_vs_alpha(all_sect_t_C, available_t_C, alpha, cL, cD, target_alpha)
    for i=1:length(all_sect_t_C)
        for j=1:size(target_alpha,2)
            try
                id(i,1) = find(all_sect_t_C(i) == available_t_C); 
                CL(i,j) = interp1(alpha,cL(:,id(i,1)),target_alpha(i,j)); % Interpolation for CL on an existing CLa curve
                CD(i,j) = interp1(alpha,cD(:,id(i,1)),target_alpha(i,j));
            catch
                if all_sect_t_C(i)>100 % Profile data correction
                    id(i,1) = find(100 == available_t_C); 
                    CL(i,j) = interp1(alpha,cL(:,id(i,1)),target_alpha(i,j)); % Interpolation for CL on 100% t/c CLa curve
                    CD(i,j) = interp1(alpha,cD(:,id(i,1)),target_alpha(i,j));
                else
                    id(i,1) = first_lower(available_t_C, all_sect_t_C(i));
                    id(i,2) = first_higher(available_t_C, all_sect_t_C(i));
                    new_CL_profi(:,1) = cL(:,id(i,1)) + (cL(:,id(i,2))-cL(:,id(i,1))).*(all_sect_t_C(i)-available_t_C(id(i,1)))./(available_t_C(id(i,2))-available_t_C(id(i,1))); % ((y-y1)/(x-x1)) = ((y2-y1)/(x2-x1))
                    new_CD_profi(:,1) = cD(:,id(i,1)) + (cD(:,id(i,2))-cD(:,id(i,1))).*(all_sect_t_C(i)-available_t_C(id(i,1)))./(available_t_C(id(i,2))-available_t_C(id(i,1)));
                    CL(i,j) = interp1(alpha,new_CL_profi,target_alpha(i,j)); % Interpolation for CL over the t/C ratio specific CLa curve
                    CD(i,j) = interp1(alpha,new_CD_profi,target_alpha(i,j));
                end
            end
        end
    end
end

function a = calc_ind_using_high_CT_approx(method, a_old, F, phi, sig, CN)
    if method == 1 % As per HANSEN eqn. 6.38
        if a_old > 0.2
            K = (4.*F.*sin(deg2rad(phi)).*sin(deg2rad(phi)))./((sig.*CN));
            a = 0.5*(2+K.*(1-2*0.2)-sqrt((((K.*(1-2*0.2))+2).^2)+4.*(K.*0.2.^2-1)));
        else
            a = 1./((4.*F.*sin(deg2rad(phi)).*sin(deg2rad(phi))./(sig.*CN))+1);
        end
    elseif method == 2 % As per HANSEN eqn 6.37
        if a_old > (1/3)
            if isnan(F)
                F=0;
            end
            p(1) = 3.*F.*sin(deg2rad(phi)).*sin(deg2rad(phi));
            p(2) = -5.*F.*sin(deg2rad(phi)).*sin(deg2rad(phi))-sig.*CN;
            p(3) = 4.*F.*sin(deg2rad(phi)).*sin(deg2rad(phi))+2.*sig.*CN;
            p(4) = -sig.*CN;
            r = roots(p);
            r=r(imag(r)==0);
            a = r(1);
        else
            a = 1./((4.*F.*sin(deg2rad(phi)).*sin(deg2rad(phi))./(sig.*CN))+1);
        end
    end
end