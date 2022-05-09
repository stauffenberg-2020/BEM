function [output_details, output, BEM] = core_bem(rho, N, op_pts, BLD , AeroFlags)
    % Core steady state BEM calculations
    % Inputs
    % 
    % rho [Kg/m^3] is the air density
    % N is the number of blades
    % wsp [m/s] is the wind speed (can be a n x 1 vector with multiple entries)
    % Pitch [deg] is the collective pitch (can be a n x 1 vector with multiple entries, But should have 1:1 correspondance with wsp)
    % rpm [rpm] is the (LSS) rotor speed in rpm (can be a nx1 vector with multiple entries, But should have 1:1 correspondance with wsp)
    % BLD is the blade data as it is from "organize_blade_data" function
    % AeroFlags are the different Aerodynamics flags

    
    % Operating data organising
    wsp = op_pts.wsp; 
    pitch = op_pts.pitch;
    rpm = op_pts.rpm;

    BEM.r = BLD.sec_r; % Blade stations
    BEM.C = BLD.sec_C; % Chord distribution
    BEM.t_C = BLD.sec_t_C; % t/C distribution
    BEM.AeroTwist = BLD.sec_twist; % Aerodynamic (mold) Twist distribution [deg]
    BEM.preflap = BLD.preflap; % Pre-bend distribution (m)
    BEM.PreBendTwist = calc_twist(BEM.r, BEM.preflap); % Pre-bend angle w.r.t root [deg]
    
    BEM.Omega_r = (rpm'*2*pi/60).*BEM.r; % Velocity seen by each of the section due to rotation [m/s]
    BEM.aA = zeros(length(BEM.r),size(wsp,1)); % Initializing axial induction factors 
    BEM.aT = zeros(length(BEM.r),size(wsp,1)); % Initializing tangential induction factors 
    BEM.V_Res = sqrt(wsp'.^2+BEM.Omega_r.^2); % Resultant velocity seen by each of the section [m/s]
    BEM.phi = atan2d((wsp'.*(1-BEM.aA)),(BEM.Omega_r.*(1+BEM.aT))); % Calculating flow angle PHI [deg]
    BEM.alpha = BEM.phi-BEM.AeroTwist-pitch'; % Angle of attack seen by each section [deg]
    [BEM.CL, BEM.CD] = CL_CD_vs_alpha(BEM.t_C, BLD.pro_t_C, BLD.pro_AoA, BLD.pro_cL, BLD.pro_cD, BEM.alpha); % CL, CD coefficients
    BEM.sig = BEM.C*N./(2*pi.*BEM.r); % Solidity ratio
    BEM.CN = BEM.CL.*cosd(BEM.phi)+BEM.CD.*sind(BEM.phi); % Sectional Normal force coefficient
    BEM.CT = BEM.CL.*sind(BEM.phi)-BEM.CD.*cosd(BEM.phi); % Sectional Tangential force coefficient
    
    if AeroFlags.induction
        for i=1:length(BEM.r) % Loop for different sections of the blade
            for j=1:size(wsp,1) % Loop for different operating set points
                BEM.aA_new(i,j) = 1./((4.*sind(BEM.phi(i,j)).*sind(BEM.phi(i,j))./(BEM.sig(i).*BEM.CN(i,j)))+1); % Calculating new axial induction factors
                BEM.aT_new(i,j) = 1./((4.*sind(BEM.phi(i,j)).*cosd(BEM.phi(i,j))./(BEM.sig(i).*BEM.CT(i,j)))-1); % Calculating new axial induction factors
                iter=1;
                while (abs(BEM.aA_new(i,j) - BEM.aA(i,j)) > 0.0001 || abs(BEM.aT_new(i,j) - BEM.aT(i,j)) > 0.0001) && iter < 1000 % Induction iteration
                    if BEM.aA_new(i,j) > 1.5, BEM.aA(i,j) = 1.5; elseif BEM.aA_new(i,j) < -1, BEM.aA(i,j) = -1; else, BEM.aA(i,j) = BEM.aA_new(i,j);end
                    if BEM.aT_new(i,j) > 1,   BEM.aT(i,j) = 1;   elseif BEM.aT_new(i,j) < -1, BEM.aT(i,j) = -1; else, BEM.aT(i,j) = BEM.aT_new(i,j);end
                    BEM.phi(i,j) = atan2d((wsp(j)'.*(1-BEM.aA(i,j))),(BEM.Omega_r(i,j).*(1+BEM.aT(i,j)))); % Calculating flow angle PHI [deg]
                    if AeroFlags.tip_loss
                        f(i,j) = (N/2).*(BEM.r(end)-BEM.r(i))./(BEM.r(i).*sind(BEM.phi(i,j)));
                        if f(i,j)<0 % Reverse flow situation, Aero twist needs to be adjusted
                            f(i,j) = 999;
                            break
                        end
                        BEM.F(i,j)=(2/pi).*acos(exp(-f(i,j)));
                    else
                        BEM.F(i,j)=1;
                    end
                    BEM.alpha(i,j) = BEM.phi(i,j)-BEM.AeroTwist(i)-pitch(j)'; % Angle of attack seen by each section [deg]
                    if BEM.alpha(i,j) < -180 || BEM.alpha(i,j) > 180
                        BEM.alpha(i,j) = mod(BEM.alpha(i,j), 180); % Check here once
                    end
                    [BEM.CL(i,j), BEM.CD(i,j)] = CL_CD_vs_alpha(BEM.t_C(i), BLD.pro_t_C, BLD.pro_AoA, BLD.pro_cL, BLD.pro_cD, BEM.alpha(i,j)); % CL, CD coefficients
                    BEM.CN(i,j) = BEM.CL(i,j).*cosd(BEM.phi(i,j))+BEM.CD(i,j).*sind(BEM.phi(i,j)); % Sectional Normal force coefficient (Thrust direction force coefficient)
                    BEM.CT(i,j) = BEM.CL(i,j).*sind(BEM.phi(i,j))-BEM.CD(i,j).*cosd(BEM.phi(i,j)); % Sectional Tangential force coefficient
                    if AeroFlags.highCT
                        method = AeroFlags.highCT; % Method = 1 or 2
                        BEM.aA_new(i,j) = calc_ind_using_high_CT_approx(method, BEM.aA(i,j), BEM.F(i,j), BEM.phi(i,j), BEM.sig(i), BEM.CN(i,j));
                    else
                        BEM.aA_new(i,j) = 1./((4.*BEM.F(i,j).*sind(BEM.phi(i,j)).*sind(BEM.phi(i,j))./(BEM.sig(i).*BEM.CN(i,j)))+1); % Calculating new axial induction factors
                    end
                    BEM.aT_new(i,j) = 1./((4.*BEM.F(i,j).*sind(BEM.phi(i,j)).*cosd(BEM.phi(i,j))./(BEM.sig(i).*BEM.CT(i,j)))-1); % Calculating new axial induction factors
                    iter=iter+1;
                end
                BEM.iter_count(i,j)=iter; % Induction iteration counter
            end
        end
    end
        
    
    %% Sectional Calculations
    A = pi.*BEM.r(end).^2; % Actual Rotor area
    Q_L = 0.5*rho.*BEM.V_Res.*BEM.V_Res; % Local Dynamic pressure for each station
    Q_G = 0.5.*rho.*wsp'.*wsp'; % Dynamic pressure
    
    % Sectional loads in Blade root coords
    BEM.FlapdF = Q_L.*BEM.C.*BEM.CN.*cosd(BEM.PreBendTwist); % Local Thrust/Normal force [N], Check this
    BEM.Rad_dF = Q_L.*BEM.C.*BEM.CN.*sind(BEM.PreBendTwist); % Local Radial force [N], Check this
    BEM.EdgedF = Q_L.*BEM.C.*BEM.CT; % Local Edge force [N]
    BEM.FlapdM = BEM.FlapdF.*(BEM.r-BEM.r(1))+BEM.Rad_dF.*(-BEM.preflap); % Local Flap Moment [Nm] at Hub end (or Blade root)
    BEM.EdgedM = BEM.EdgedF.*(BEM.r-BEM.r(1)); % Local Edge Moment [Nm] at Hub end (or Blade root)
    
    % Sectional loads w.r.t rotor center
    BEM.EdgedMC = BEM.EdgedF.*BEM.r; % Local Edge Moment [Nm] (w.r.t rotor center) in blade root coords
    
    % Full Rotor loads and coefficients
    BEM.dT = 3.*BEM.FlapdF;
    BEM.Ct = BEM.dT./(Q_G.*2.*pi.*BEM.r); % Local Thrust Coefficient from local Thrust
    

    
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
            K = (4.*F.*sind(phi).*sind(phi))./((sig.*CN));
            a = 0.5*(2+K.*(1-2*0.2)-sqrt((((K.*(1-2*0.2))+2).^2)+4.*(K.*0.2.^2-1)));
        else
            a = 1./((4.*F.*sind(phi).*sind(phi)./(sig.*CN))+1);
        end
    elseif method == 2 % As per HANSEN eqn 6.37
        if a_old > (1/3)
            if isnan(F)
                F=0;
            end
            p(1) = 3.*F.*sind(phi).*sind(phi);
            p(2) = -5.*F.*sind(phi).*sind(phi)-sig.*CN;
            p(3) = 4.*F.*sind(phi).*sind(phi)+2.*sig.*CN;
            p(4) = -sig.*CN;
            r = roots(p);
            r=r(imag(r)==0);
            [~,idx] = min(abs(r-a_old)); % Finding the positive root near to previous 'a'
            a = r(idx);
        else
            a = 1./((4.*F.*sind(phi).*sind(phi)./(sig.*CN))+1);
        end
    end
end

function bend_angle = calc_twist(sec_R, position_R)
    bend_angle = atand(-position_R./sec_R); % Calculating the pre bend angle (relative to root) due to given position, Is this accurate/correct?
end