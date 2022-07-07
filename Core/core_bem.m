 function [output_details, output, BEM] = core_bem(General, op_pts, BLD)
    % Core steady state BEM calculations
    % Comments
    % 
    % 
    
    % Operating data organising
    wsp = op_pts.wsp; 
    pitch = op_pts.pitch;
    rpm = op_pts.rpm;
    
    % Pre-initializing structure for MATLAB Codegen
    BEM = struct('r',0,'C',0,'t_C',0,'AeroTwist',0,'preflap',0,'PreBendTwist',0,...
                 'Omega_r',0,'aA',0,'aA_new',0,'aT',0,'aT_new',0,'V_Res',0,'phi',0,'alpha',0,'CL',0,'CD',0,...
                 'sig',0,'CN',0,'CT',0,'f',0,'F',0,'iter_count',0,...
                 'FlapdF',0,'Rad_dF',0,'EdgedF',0,'FlapdM',0,'EdgedM',0,'EdgedMC',0,'dT',0,'Ct',0);
    coder.varsize('BEM.r','BEM.C','BEM.t_C','BEM.AeroTwist','BEM.preflap','BEM.PreBendTwist',...
                  'BEM.Omega_r','BEM.aA','BEM.aA_new','BEM.aT','BEM.aT_new','BEM.V_Res','BEM.phi','BEM.alpha','BEM.CL','BEM.CD',...
                  'BEM.sig','BEM.CN','BEM.CT','BEM.f','BEM.F','BEM.iter_count',...
                  'BEM.FlapdF','BEM.Rad_dF','BEM.EdgedF','BEM.FlapdM','BEM.EdgedM','BEM.EdgedMC','BEM.dT','BEM.Ct'); 
    
    BEM.r = BLD.r; % Blade stations
    BEM.C = BLD.C; % Chord distribution
    BEM.t_C = BLD.t_C; % t/C distribution
    BEM.AeroTwist = BLD.AeroTwist; % Aerodynamic (mold) Twist distribution [deg]
    BEM.preflap = BLD.preflap; % Pre-bend distribution (m)
    BEM.PreBendTwist = calc_twist(BEM.r, BEM.preflap); % Pre-bend angle w.r.t root [deg]
    
    BEM.Omega_r = (repmat(rpm',length(BEM.r),1)*2*pi/60).*repmat(BEM.r,1,length(rpm)); % Velocity seen by each of the section due to rotation [m/s]
    BEM.aA = zeros(length(BEM.r),size(wsp,1)); % Initializing axial induction factors 
    BEM.aT = zeros(length(BEM.r),size(wsp,1)); % Initializing tangential induction factors 
    BEM.V_Res = sqrt(repmat((wsp.^2)',length(BEM.r),1)+BEM.Omega_r.^2); % Resultant velocity seen by each of the section [m/s]
    BEM.phi = atan2d((repmat(wsp',length(BEM.r),1).*(1-BEM.aA)),(BEM.Omega_r.*(1+BEM.aT))); % Calculating flow angle PHI [deg]
    BEM.alpha = BEM.phi-repmat(BEM.AeroTwist,1,length(wsp))-repmat(pitch',length(BEM.r),1); % Angle of attack seen by each section [deg]
    [BEM.CL, BEM.CD] = CL_CD_vs_alpha(BEM.t_C, BLD.pro_t_C, BLD.pro_AoA, BLD.pro_cL, BLD.pro_cD, BEM.alpha); % CL, CD coefficients
    BEM.sig = BEM.C*General.N./(2*pi.*BEM.r); % Solidity ratio
    BEM.CN = BEM.CL.*cosd(BEM.phi)+BEM.CD.*sind(BEM.phi); % Sectional Normal force coefficient
    BEM.CT = BEM.CL.*sind(BEM.phi)-BEM.CD.*cosd(BEM.phi); % Sectional Tangential force coefficient
    
    if General.induction
        BEM.aA_new = zeros(length(BEM.r),size(wsp,1)); % Initializing axial induction factors 
        BEM.aT_new = zeros(length(BEM.r),size(wsp,1)); % Initializing tangential induction factors 
        BEM.f = zeros(length(BEM.r),size(wsp,1)); % Initializing f matrix
        BEM.F = zeros(length(BEM.r),size(wsp,1)); % Initializing f matrix
        BEM.iter_count = zeros(length(BEM.r),size(wsp,1)); % Initializing f matrix
        for i=1:length(BEM.r) % Loop for different sections of the blade
            for j=1:size(wsp,1) % Loop for different operating set points
                BEM.aA_new(i,j) = 1./((4.*sind(BEM.phi(i,j)).*sind(BEM.phi(i,j))./(BEM.sig(i).*BEM.CN(i,j)))+1); % Calculating new axial induction factors
                BEM.aT_new(i,j) = 1./((4.*sind(BEM.phi(i,j)).*cosd(BEM.phi(i,j))./(BEM.sig(i).*BEM.CT(i,j)))-1); % Calculating new axial induction factors
                iter=1;
                while (abs(BEM.aA_new(i,j) - BEM.aA(i,j)) > 0.0001 || abs(BEM.aT_new(i,j) - BEM.aT(i,j)) > 0.0001) && iter < 1000 % Induction iteration
                    if BEM.aA_new(i,j) > 1.5, BEM.aA(i,j) = 1.5; elseif BEM.aA_new(i,j) < -1, BEM.aA(i,j) = -1; else, BEM.aA(i,j) = BEM.aA_new(i,j);end
                    if BEM.aT_new(i,j) > 1,   BEM.aT(i,j) = 1;   elseif BEM.aT_new(i,j) < -1, BEM.aT(i,j) = -1; else, BEM.aT(i,j) = BEM.aT_new(i,j);end
                    BEM.phi(i,j) = atan2d((wsp(j)'.*(1-BEM.aA(i,j))),(BEM.Omega_r(i,j).*(1+BEM.aT(i,j)))); % Calculating flow angle PHI [deg]
                    if General.tip_loss
                        BEM.f(i,j) = (General.N/2).*(BEM.r(end)-BEM.r(i))./(BEM.r(i).*sind(BEM.phi(i,j)));
                        if BEM.f(i,j)<0 % Reverse flow situation, Aero twist needs to be adjusted
                            BEM.f(i,j) = 999;
                        end
                        BEM.F(i,j)=(2/pi).*acos(exp(-BEM.f(i,j)));
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
                    if General.highCT
                        method = General.highCT; % Method = 1 or 2
                        BEM.aA_new(i,j) = real(calc_ind_using_high_CT_approx(method, BEM.aA(i,j), BEM.F(i,j), BEM.phi(i,j), BEM.sig(i), BEM.CN(i,j)));
                    else
                        BEM.aA_new(i,j) = 1./((4.*BEM.F(i,j).*sind(BEM.phi(i,j)).*sind(BEM.phi(i,j))./(BEM.sig(i).*BEM.CN(i,j)))+1); % Calculating new axial induction factors
                    end
                    BEM.aT_new(i,j) = 1./((4.*BEM.F(i,j).*sind(BEM.phi(i,j)).*cosd(BEM.phi(i,j))./(BEM.sig(i).*BEM.CT(i,j)))-1); % Calculating new axial induction factors
                    RF=0.50; % Relaxation Factor as per https://onlinelibrary.wiley.com/doi/full/10.1002/ese3.945
                    BEM.aA_new(i,j) = RF*BEM.aA_new(i,j)+(1-RF)*BEM.aA(i,j);
                    BEM.aT_new(i,j) = RF*BEM.aT_new(i,j)+(1-RF)*BEM.aT(i,j);
                    iter=iter+1;
                end
                BEM.iter_count(i,j)=iter; % Induction iteration counter
            end
        end
    end
        
    
    %% Sectional Calculations
    A = pi.*BEM.r(end).^2; % Actual Rotor area
    Q_L = 0.5*General.rho.*BEM.V_Res.*BEM.V_Res; % Local Dynamic pressure for each station
    Q_G = 0.5.*General.rho.*wsp'.*wsp'; % Dynamic pressure
    
    % Sectional loads in Blade root coords
    BEM.FlapdF = Q_L.*repmat(BEM.C,1,length(wsp)).*BEM.CN.*repmat(cosd(BEM.PreBendTwist),1,length(wsp)); % Local Thrust/Normal force [N], Check this
    BEM.Rad_dF = Q_L.*repmat(BEM.C,1,length(wsp)).*BEM.CN.*repmat(sind(BEM.PreBendTwist),1,length(wsp)); % Local Radial force [N], Check this
    BEM.EdgedF = Q_L.*repmat(BEM.C,1,length(wsp)).*BEM.CT; % Local Edge force [N]
    BEM.FlapdM = BEM.FlapdF.*repmat(BEM.r-BEM.r(1),1,length(wsp))+BEM.Rad_dF.*repmat((-BEM.preflap),1,length(wsp)); % Local Flap Moment [Nm] at Hub end (or Blade root)
    BEM.EdgedM = BEM.EdgedF.*repmat(BEM.r-BEM.r(1),1,length(wsp)); % Local Edge Moment [Nm] at Hub end (or Blade root)
    
    % Sectional loads w.r.t rotor center
    BEM.EdgedMC = BEM.EdgedF.*repmat(BEM.r,1,length(wsp)); % Local Edge Moment [Nm] (w.r.t rotor center) in blade root coords
    
    % Full Rotor loads and coefficients
    BEM.dT = General.N.*BEM.FlapdF;
    BEM.Ct = BEM.dT./(repmat(Q_G,length(BEM.r),1).*2.*pi.*repmat(BEM.r,1,length(wsp))); % Local Thrust Coefficient from local Thrust
    
    %% Output calculations
    output_details = {'Flp_F', 'Edg_F', 'Flp_M', 'Edg_M', 'M_aer', 'P_rot', 'Cp', 'Fthr', 'Ct';
                            'Flap_For', 'Edge_For', 'Flap_Mom', 'Edge_Mom', 'Aero_Tor', 'Aero_Pow', 'Rotor_Cp','Aero_Thr','Rotor_Ct';
                            '[kN]','[kN]','[kNm]','[kNm]','[kNm]','[kW]','[-]','[kN]','[-]'};
    
    output=zeros(length(wsp),length(output_details));
    
    output(1:length(wsp),1) = (1/1000)*trapz(BEM.r(1:length(BEM.r)-1),BEM.FlapdF(1:length(BEM.r)-1,1:length(wsp))); % Flap Force [KN], % Ignoring the influence of last section similar to VTS
    output(1:length(wsp),2) = (1/1000)*trapz(BEM.r(1:length(BEM.r)-1),BEM.EdgedF(1:length(BEM.r)-1,1:length(wsp))); % Edge Force [KN]
    output(1:length(wsp),3) = (1/1000)*trapz(BEM.r(1:length(BEM.r)-1),BEM.FlapdM(1:length(BEM.r)-1,1:length(wsp))); % Flap Moment [KNm]
    output(1:length(wsp),4) = (1/1000)*trapz(BEM.r(1:length(BEM.r)-1),BEM.EdgedM(1:length(BEM.r)-1,1:length(wsp))); % Edge Moment [KNm]
    output(1:length(wsp),5) = (1/1000)*trapz(BEM.r(1:length(BEM.r)-1),BEM.EdgedMC(1:length(BEM.r)-1,1:length(wsp)))*General.N; % Aerodynamic Torque [KNm]
    output(1:length(wsp),6) = output(:,5).*(2*pi.*rpm/60); % Rotor Power [KW]
    output(1:length(wsp),7) = (output(:,6).*1000)./(Q_G'.*A.*wsp); % Rotor averaged Cp calculated from rotor power
    output(1:length(wsp),8) = output(:,1)*General.N; % Aerodynamic Thrust [KN]
    output(1:length(wsp),9) = (output(:,8).*1000)./(Q_G'.*A); % Rotor averaged Ct Calculated from rotor thrust
    
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
    id = zeros(length(all_sect_t_C),2);
    CL = zeros(length(all_sect_t_C),size(target_alpha,2));
    CD = zeros(length(all_sect_t_C),size(target_alpha,2));
    new_CL_profi = zeros(size(cL,1),1);
    new_CD_profi = zeros(size(cL,1),1);
    for i=1:length(all_sect_t_C)
        for j=1:size(target_alpha,2)
            if any(all_sect_t_C(i) == available_t_C)
                index = find(all_sect_t_C(i) == available_t_C); 
                id(i,1) = index(1); % Multi line statement for coder compatibility
                CL(i,j) = interp1(alpha,cL(:,id(i,1)),target_alpha(i,j)); % Interpolation for CL on an existing CLa curve
                CD(i,j) = interp1(alpha,cD(:,id(i,1)),target_alpha(i,j));
            elseif all_sect_t_C(i)>100 % Profile data correction
                index = find(100 == available_t_C);
                id(i,1) = index(1); % Multi line statement for coder compatibility
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

function a = calc_ind_using_high_CT_approx(method, a_old, F, phi, sig, CN)
    a = 0; % Initializing variable for coder compatibility
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
            p=zeros(1,4);
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