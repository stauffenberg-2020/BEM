 function [output_details, output, BEM] = core_bem(General, op_pts, BLD)
    % Core steady state BEM calculations
    % Comments
    % 
    % 
    
    % Operating data organising
    wsp = op_pts.wsp;
    shear = op_pts.shear;
    pitch = op_pts.pitch;
    rpm = op_pts.rpm;
    
    if all(shear == 0), bld_calcs = 1; else bld_calcs =3; end
    
%     BEM = initialize_structs(); % Pre-initializing structure for MATLAB Codegen
              
    BEM.r = BLD.r; % Blade stations
    BEM.C = BLD.C; % Chord distribution
    BEM.t_C = BLD.t_C; % t/C distribution
    BEM.AeroTwist = BLD.AeroTwist; % Aerodynamic (mold) Twist distribution [deg]
    BEM.preflap = BLD.preflap; % Pre-bend distribution (m)
    PRO.AoA = BLD.pro_AoA;
    PRO.t_C = BLD.pro_t_C;
    PRO.CL = BLD.pro_cL;
    PRO.CD = BLD.pro_cD;
    
    BEM.PreBendTwist = calc_twist(BEM.r, BEM.preflap); % Pre-bend angle w.r.t root [deg]
    BEM.Omega_r = (repmat(rpm',length(BEM.r),1)*2*pi/60).*repmat(BEM.r,1,length(rpm)); % Velocity seen by each of the section due to rotation [m/s]
    BEM.sig = BEM.C*General.N./(2*pi.*BEM.r); % Solidity ratio
    
    bld_id = {'b1','b2','b3'}; % b = fieldnames(BEM.aA); % Getting the blade struct field names
    if bld_calcs == 1 % Condition for whether to skip 2nd and 3rd blade calculation
        id = 1;
        b = bld_id{id};
        az = ['az_ind']; 
        BEM.(b).(az).aA = zeros(length(BEM.r),size(wsp,1)); % Initializing axial induction factors 
        BEM.(b).(az).aT = zeros(length(BEM.r),size(wsp,1)); % Initializing tangential induction factors 
        BEM.(b).(az).V_Res = sqrt(repmat((wsp.^2)',length(BEM.r),1)+BEM.Omega_r.^2); % Resultant velocity seen by each of the section [m/s]
        BEM.(b).(az).phi = atan2d((repmat(wsp',length(BEM.r),1).*(1-BEM.(b).(az).aA)),(BEM.Omega_r.*(1+BEM.(b).(az).aT))); % Calculating flow angle PHI [deg]
        BEM.(b).(az).alpha = BEM.(b).(az).phi-repmat(BEM.AeroTwist,1,length(wsp))-repmat(pitch',length(BEM.r),1); % Angle of attack seen by each section [deg]
        [BEM.(b).(az).CL, BEM.(b).(az).CD] = CL_CD_vs_alpha(BEM.t_C, BLD.pro_t_C, BLD.pro_AoA, BLD.pro_cL, BLD.pro_cD, BEM.(b).(az).alpha); % CL, CD coefficients
        BEM.(b).(az).CN = BEM.(b).(az).CL.*cosd(BEM.(b).(az).phi)+BEM.(b).(az).CD.*sind(BEM.(b).(az).phi); % Sectional Normal force coefficient
        BEM.(b).(az).CT = BEM.(b).(az).CL.*sind(BEM.(b).(az).phi)-BEM.(b).(az).CD.*cosd(BEM.(b).(az).phi); % Sectional Tangential force coefficient
        if General.induction % Condition for doing induction calculations
            tmp_BEM_inp = BEM_to_temp_BEM(BEM, b, az, 1:length(BEM.r), 1:length(wsp));
            BLD.R = BLD.r(end);
            BEM.(b).(az) = iterate_for_induction(General, tmp_BEM_inp, PRO, BLD, wsp, pitch);
            clear tmp_BEM_inp;
        end

    else
        for id=1:numel(bld_id) % Loop for each of the blades
            b = bld_id{id};
            azimuths = [0:30:360];
            for k = 1:length(azimuths) % Loop for different azimuth positions of blade
                az = ['az_',num2str(azimuths(k),'%03d')];
                BEM.(b).(az).aA = zeros(length(BEM.r),size(wsp,1)); % Initializing axial induction factors 
                BEM.(b).(az).aT = zeros(length(BEM.r),size(wsp,1)); % Initializing tangential induction factors 
                for o = 1:length(wsp) % Loop for different operating points
                    for s=1:length(BEM.r) % Loop for different sections of the blade
                        if id == 2, azi = mod(azimuths(k)+120,360); elseif id ==3, azi = mod(azimuths(k)+240,360); else azi = mod(azimuths(k),360); end;
                        BEM.(b).(az).wsp_sh(s,o) = wsp_vs_azim(wsp(o), shear(o), azi ,90, BEM.r(s)); % sheared wind speed across height w.r.t rotor azimuth
                        BEM.(b).(az).V_Res(s,o) = sqrt(repmat((BEM.(b).(az).wsp_sh(s,o).^2)',length(BEM.r(s)),1)+BEM.Omega_r(s,o).^2); % Resultant velocity seen by each of the section [m/s]
                        BEM.(b).(az).phi(s,o) = atan2d((repmat(BEM.(b).(az).wsp_sh(s,o)',length(BEM.r(s)),1).*(1-BEM.(b).(az).aA(s,o))),(BEM.Omega_r(s,o).*(1+BEM.(b).(az).aT(s,o)))); % Calculating flow angle PHI [deg]
                        BEM.(b).(az).alpha(s,o) = BEM.(b).(az).phi(s,o)-BEM.AeroTwist(s)-pitch(o); % Angle of attack seen by each section [deg]
                        [BEM.(b).(az).CL(s,o), BEM.(b).(az).CD(s,o)] = CL_CD_vs_alpha(BEM.t_C(s), BLD.pro_t_C, BLD.pro_AoA, BLD.pro_cL, BLD.pro_cD, BEM.(b).(az).alpha(s,o)); % CL, CD coefficients
                        BEM.(b).(az).CN(s,o) = BEM.(b).(az).CL(s,o).*cosd(BEM.(b).(az).phi(s,o))+BEM.(b).(az).CD(s,o).*sind(BEM.(b).(az).phi(s,o)); % Sectional Normal force coefficient
                        BEM.(b).(az).CT(s,o) = BEM.(b).(az).CL(s,o).*sind(BEM.(b).(az).phi(s,o))-BEM.(b).(az).CD(s,o).*cosd(BEM.(b).(az).phi(s,o)); % Sectional Tangential force coefficient
                        if General.induction % Condition for doing induction calculations
                            tmp_BEM_inp = BEM_to_temp_BEM(BEM, b, az, s, o);
                            local_BLD = BLD_to_local_BLD(BLD, s);
                            tmp_BEM_out = iterate_for_induction(General, tmp_BEM_inp, PRO, local_BLD, BEM.(b).(az).wsp_sh(s,o), pitch(o));
                            BEM = temp_BEM_to_BEM(tmp_BEM_out, BEM, b, az, s, o);
                            clear tmp_inp_BEM local_BLD tmp_BEM_out;
                        end
                    end
                end
            end
        end
    end
    
    BEM = average_across_azimuth(BEM);
        
    
    %% Sectional Calculations
    A = pi.*BEM.r(end).^2; % Actual Rotor area
    Q_G = 0.5.*General.rho.*wsp'.*wsp'; % Dynamic pressure
    
    % Sectional loads in Blade root coords
    for id=1:numel(bld_id)
        b = bld_id{id};
        if id <= bld_calcs
            BEM.(b).Q_L = 0.5*General.rho.*BEM.(b).V_Res_avg.*BEM.(b).V_Res_avg; % Local Dynamic pressure for each station
            BEM.(b).FlapdF = BEM.(b).Q_L.*repmat(BEM.C,1,length(wsp)).*BEM.(b).CN_avg.*repmat(cosd(BEM.PreBendTwist),1,length(wsp)); % Local Thrust/Normal force [N], Check this
            BEM.(b).Rad_dF = BEM.(b).Q_L.*repmat(BEM.C,1,length(wsp)).*BEM.(b).CN_avg.*repmat(sind(BEM.PreBendTwist),1,length(wsp)); % Local Radial force [N], Check this
            BEM.(b).EdgedF = BEM.(b).Q_L.*repmat(BEM.C,1,length(wsp)).*BEM.(b).CT_avg; % Local Edge force [N]
            BEM.(b).FlapdM = BEM.(b).FlapdF.*repmat(BEM.r-BEM.r(1),1,length(wsp))+BEM.(b).Rad_dF.*repmat((-BEM.preflap),1,length(wsp)); % Local Flap Moment [Nm] at Hub end (or Blade root)
            BEM.(b).EdgedM = BEM.(b).EdgedF.*repmat(BEM.r-BEM.r(1),1,length(wsp)); % Local Edge Moment [Nm] at Hub end (or Blade root)
        end
    end
    
    if bld_calcs == 1
        BEM.FlapdF_max = BEM.b1.FlapdF;
        BEM.EdgedF_max = BEM.b1.EdgedF;
        BEM.FlapdM_max = BEM.b1.FlapdM;
        BEM.EdgedM_max = BEM.b1.EdgedM;
        BEM.EdgedMC_all = General.N*BEM.b1.EdgedF.*repmat(BEM.r,1,length(wsp));
        BEM.dT = General.N*BEM.b1.FlapdF;
    else
        BEM.FlapdF_max = max(max(BEM.b1.FlapdF, BEM.b2.FlapdF), BEM.b3.FlapdF);
        BEM.EdgedF_max = max(max(BEM.b1.EdgedF, BEM.b2.EdgedF), BEM.b3.EdgedF);
        BEM.FlapdM_max = max(max(BEM.b1.FlapdM, BEM.b2.FlapdM), BEM.b3.FlapdM);
        BEM.EdgedM_max = max(max(BEM.b1.EdgedM, BEM.b2.EdgedM), BEM.b3.EdgedM);
        BEM.EdgedMC_all = (BEM.b1.EdgedF + BEM.b2.EdgedF + BEM.b3.EdgedF).*repmat(BEM.r,1,length(wsp)); % Local Edge Moment [Nm] (w.r.t rotor center) in blade root coords
        BEM.dT = BEM.b1.FlapdF+BEM.b2.FlapdF+BEM.b3.FlapdF; % Full Rotor loads
    end 
    BEM.Ct = BEM.dT./(repmat(Q_G,length(BEM.r),1).*2.*pi.*repmat(BEM.r,1,length(wsp))); % Local Thrust Coefficient from local Thrust
    
    %% Output calculations
    output_details = {'Flp_F', 'Edg_F', 'Flp_M', 'Edg_M', 'M_aer', 'P_rot', 'Cp', 'Fthr', 'Ct';
                      'Flap_For', 'Edge_For', 'Flap_Mom', 'Edge_Mom', 'Aero_Tor', 'Aero_Pow', 'Rotor_Cp','Aero_Thr','Rotor_Ct';
                      '[kN]','[kN]','[kNm]','[kNm]','[kNm]','[kW]','[-]','[kN]','[-]'};
    
    output=zeros(length(wsp),length(output_details));
    
    output(1:length(wsp),1) = (1/1000)*trapz(BEM.r(1:length(BEM.r)-1),BEM.FlapdF_max(1:length(BEM.r)-1,1:length(wsp))); % Flap Force [KN], % Ignoring the influence of last section similar to VTS
    output(1:length(wsp),2) = (1/1000)*trapz(BEM.r(1:length(BEM.r)-1),BEM.EdgedF_max(1:length(BEM.r)-1,1:length(wsp))); % Edge Force [KN]
    output(1:length(wsp),3) = (1/1000)*trapz(BEM.r(1:length(BEM.r)-1),BEM.FlapdM_max(1:length(BEM.r)-1,1:length(wsp))); % Flap Moment [KNm]
    output(1:length(wsp),4) = (1/1000)*trapz(BEM.r(1:length(BEM.r)-1),BEM.EdgedM_max(1:length(BEM.r)-1,1:length(wsp))); % Edge Moment [KNm]
    output(1:length(wsp),5) = (1/1000)*trapz(BEM.r(1:length(BEM.r)-1),BEM.EdgedMC_all(1:length(BEM.r)-1,1:length(wsp))); % Aerodynamic Torque [KNm]
    output(1:length(wsp),6) = output(:,5).*(2*pi.*rpm/60); % Rotor Power [KW]
    output(1:length(wsp),7) = (output(:,6).*1000)./(Q_G'.*A.*wsp); % Rotor averaged Cp calculated from rotor power
    output(1:length(wsp),8) = (1/1000)*trapz(BEM.r(1:length(BEM.r)-1),BEM.dT(1:length(BEM.r)-1,1:length(wsp))); % Aerodynamic Thrust [KN]
    output(1:length(wsp),9) = (output(:,8).*1000)./(Q_G'.*A); % Rotor averaged Ct Calculated from rotor thrust
    
end

%% Supporting functions
function BEM = iterate_for_induction(General, BEM, PRO, BLD, wsp, pitch)
    BEM.aA_new = zeros(length(BLD.r),size(wsp,1)); % Initializing axial induction factors 
    BEM.aT_new = zeros(length(BLD.r),size(wsp,1)); % Initializing tangential induction factors 
    BEM.f = zeros(length(BLD.r),size(wsp,1)); % Initializing f matrix
    BEM.F = zeros(length(BLD.r),size(wsp,1)); % Initializing F matrix
    BEM.iter_count = zeros(length(BLD.r),size(wsp,1)); % Initializing iter count matrix
    for s=1:length(BLD.r) % Loop for different sections of the blade
        for o=1:size(wsp,1) % Loop for different Operating points
            BEM.aA_new(s,o) = 1./((4.*sind(BEM.phi(s,o)).*sind(BEM.phi(s,o))./(BEM.sig(s).*BEM.CN(s,o)))+1); % Calculating new axial induction factors
            BEM.aT_new(s,o) = 1./((4.*sind(BEM.phi(s,o)).*cosd(BEM.phi(s,o))./(BEM.sig(s).*BEM.CT(s,o)))-1); % Calculating new axial induction factors
            iter=1;
            while (abs(BEM.aA_new(s,o) - BEM.aA(s,o)) > 0.0001 || abs(BEM.aT_new(s,o) - BEM.aT(s,o)) > 0.0001) && iter < 1000 % Induction iteration
                if BEM.aA_new(s,o) > 1.5, BEM.aA(s,o) = 1.5; elseif BEM.aA_new(s,o) < -1, BEM.aA(s,o) = -1; else, BEM.aA(s,o) = BEM.aA_new(s,o);end
                if BEM.aT_new(s,o) > 1,   BEM.aT(s,o) = 1;   elseif BEM.aT_new(s,o) < -1, BEM.aT(s,o) = -1; else, BEM.aT(s,o) = BEM.aT_new(s,o);end
                BEM.phi(s,o) = atan2d((wsp(o)'.*(1-BEM.aA(s,o))),(BEM.Omega_r(s,o).*(1+BEM.aT(s,o)))); % Calculating flow angle PHI [deg]
                if General.tip_loss
                    BEM.f(s,o) = (General.N/2).*(BLD.R-BLD.r(s))./(BLD.r(s).*sind(BEM.phi(s,o)));
                    if BEM.f(s,o)<0 % Reverse flow situation, Aero twist needs to be adjusted
                        BEM.f(s,o) = 999;
                    end
                    BEM.F(s,o)=(2/pi).*acos(exp(-BEM.f(s,o)));
                else
                    BEM.F(s,o)=1;
                end
                BEM.alpha(s,o) = BEM.phi(s,o)-BLD.AeroTwist(s)-pitch(o)'; % Angle of attack seen by each section [deg]
                if BEM.alpha(s,o) < -180 || BEM.alpha(s,o) > 180
                    BEM.alpha(s,o) = mod(BEM.alpha(s,o), 180); % Check here once
                end
                [BEM.CL(s,o), BEM.CD(s,o)] = CL_CD_vs_alpha(BLD.t_C(s), PRO.t_C, PRO.AoA, PRO.CL, PRO.CD, BEM.alpha(s,o)); % CL, CD coefficients
                BEM.CN(s,o) = BEM.CL(s,o).*cosd(BEM.phi(s,o))+BEM.CD(s,o).*sind(BEM.phi(s,o)); % Sectional Normal force coefficient (Thrust direction force coefficient)
                BEM.CT(s,o) = BEM.CL(s,o).*sind(BEM.phi(s,o))-BEM.CD(s,o).*cosd(BEM.phi(s,o)); % Sectional Tangential force coefficient
                if General.highCT
                    method = General.highCT; % Method = 1 or 2
                    BEM.aA_new(s,o) = real(calc_ind_using_high_CT_approx(method, BEM.aA(s,o), BEM.F(s,o), BEM.phi(s,o), BEM.sig(s), BEM.CN(s,o)));
                else
                    BEM.aA_new(s,o) = 1./((4.*BEM.F(s,o).*sind(BEM.phi(s,o)).*sind(BEM.phi(s,o))./(BEM.sig(s).*BEM.CN(s,o)))+1); % Calculating new axial induction factors
                end
                BEM.aT_new(s,o) = 1./((4.*BEM.F(s,o).*sind(BEM.phi(s,o)).*cosd(BEM.phi(s,o))./(BEM.sig(s).*BEM.CT(s,o)))-1); % Calculating new axial induction factors
                RF=0.50; % Relaxation Factor as per https://onlinelibrary.wiley.com/doi/full/10.1002/ese3.945
                BEM.aA_new(s,o) = RF*BEM.aA_new(s,o)+(1-RF)*BEM.aA(s,o);
                BEM.aT_new(s,o) = RF*BEM.aT_new(s,o)+(1-RF)*BEM.aT(s,o);
                iter=iter+1;
            end
            BEM.iter_count(s,o)=iter; % Induction iteration counter
        end
    end
end

function BEM = average_across_azimuth(BEM)
    
    az = fieldnames(BEM.b1);
    if length(az) == 1
        BEM.b1.CN_avg = BEM.b1.(az{1}).CN;
        BEM.b1.CT_avg = BEM.b1.(az{1}).CT;
        BEM.b1.V_Res_avg = BEM.b1.(az{1}).V_Res;
    else
        % Averaging across rotor azimuths
        bld_id = {'b1','b2','b3'};
        for id =1:numel(bld_id)
            b = bld_id{id};
            BEM.(b).CN_avg = zeros(size(BEM.Omega_r));
            BEM.(b).CT_avg = zeros(size(BEM.Omega_r));
            BEM.(b).V_Res_avg = zeros(size(BEM.Omega_r));
            for k = 1:length(az)
                
                BEM.(b).CN_avg = BEM.(b).CN_avg + BEM.(b).(az{k}).CN;
                BEM.(b).CT_avg = BEM.(b).CT_avg + BEM.(b).(az{k}).CT;
                BEM.(b).V_Res_avg = BEM.(b).V_Res_avg + BEM.(b).(az{k}).V_Res;
            end
        end

        for id = 1:numel(bld_id)
            b = bld_id{id};
            BEM.(b).CN_avg = BEM.(b).CN_avg/length(az);
            BEM.(b).CT_avg = BEM.(b).CT_avg/length(az);
            BEM.(b).V_Res_avg = BEM.(b).V_Res_avg/length(az);
        end
    end
    
end

function tmp_BEM = BEM_to_temp_BEM(BEM, b, az, s, o)
    tmp_BEM.aA = BEM.(b).(az).aA(s,o);
    tmp_BEM.aT = BEM.(b).(az).aT(s,o);
    tmp_BEM.V_Res = BEM.(b).(az).V_Res(s,o);
    tmp_BEM.phi = BEM.(b).(az).phi(s,o);
    tmp_BEM.alpha = BEM.(b).(az).alpha(s,o);
    tmp_BEM.CL = BEM.(b).(az).CL(s,o);
    tmp_BEM.CD = BEM.(b).(az).CD(s,o);
    tmp_BEM.CN = BEM.(b).(az).CN(s,o);
    tmp_BEM.CT = BEM.(b).(az).CT(s,o);
    tmp_BEM.Omega_r = BEM.Omega_r(s,o);
    tmp_BEM.sig = BEM.sig(s);
end

function BEM = temp_BEM_to_BEM(tmp_BEM, BEM, b, az, s, o)
    BEM.(b).(az).aA(s,o) = tmp_BEM.aA;
    BEM.(b).(az).aT(s,o) = tmp_BEM.aT;
    BEM.(b).(az).V_Res(s,o) = tmp_BEM.V_Res;
    BEM.(b).(az).phi(s,o) = tmp_BEM.phi;
    BEM.(b).(az).alpha(s,o) = tmp_BEM.alpha;
    BEM.(b).(az).CL(s,o) = tmp_BEM.CL;
    BEM.(b).(az).CD(s,o) = tmp_BEM.CD;
    BEM.(b).(az).CN(s,o) = tmp_BEM.CN;
    BEM.(b).(az).CT(s,o) = tmp_BEM.CT;
end

function local_BLD = BLD_to_local_BLD(BLD, s)
    local_BLD.r = BLD.r(s);
    local_BLD.AeroTwist = BLD.AeroTwist(s);
    local_BLD.t_C = BLD.t_C(s);
    local_BLD.R = BLD.r(end);
end

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

function sheared_wsp = wsp_vs_azim(wsp, alpha, azim, ref_height, bld_length)
    heights = ref_height-bld_length*(cosd(azim));
    sheared_wsp = wsp.*(heights./ref_height).^(alpha);
end

function bend_angle = calc_twist(sec_R, position_R)
    bend_angle = atand(-position_R./sec_R); % Calculating the pre bend angle (relative to root) due to given position, Is this accurate/correct?
end

function BEM = initialize_structs()
    BEM = struct('r',0,'C',0,'t_C',0,'AeroTwist',0,'preflap',0,'PreBendTwist',0,'Omega_r',0,...
             'aA',struct('b1',0,'b2',0,'b3',0),'aA_new',struct('b1',0,'b2',0,'b3',0),...
             'aT',struct('b1',0,'b2',0,'b3',0),'aT_new',struct('b1',0,'b2',0,'b3',0),...
             'V_Res',struct('b1',0,'b2',0,'b3',0),'phi',struct('b1',0,'b2',0,'b3',0),...
             'alpha',struct('b1',0,'b2',0,'b3',0),'CL',struct('b1',0,'b2',0,'b3',0),...
             'CD',struct('b1',0,'b2',0,'b3',0),'sig',0,...
             'CN',struct('b1',0,'b2',0,'b3',0),'CT',struct('b1',0,'b2',0,'b3',0),...
             'f',struct('b1',0,'b2',0,'b3',0),'F',struct('b1',0,'b2',0,'b3',0),...
             'iter_count',struct('b1',0,'b2',0,'b3',0),'FlapdF',struct('b1',0,'b2',0,'b3',0),...
             'Rad_dF',struct('b1',0,'b2',0,'b3',0),'EdgedF',struct('b1',0,'b2',0,'b3',0),...
             'FlapdM',struct('b1',0,'b2',0,'b3',0),'EdgedM',struct('b1',0,'b2',0,'b3',0),...
             'EdgedMC',struct('b1',0,'b2',0,'b3',0),'dT',0,'Ct',0,'Q_L',struct('b1',0,'b2',0,'b3',0),...
             'FlapdF_max',0,'EdgedF_max',0,'FlapdM_max',0,'EdgedM_max',0,'EdgedMC_all',0);
    coder.varsize('BEM.r','BEM.C','BEM.t_C','BEM.AeroTwist','BEM.preflap','BEM.PreBendTwist','BEM.Omega_r',...
                  'BEM.aA.b1','BEM.aA.b2','BEM.aA.b3','BEM.aA_new.b1','BEM.aA_new.b2','BEM.aA_new.b3',...
                  'BEM.aT.b1','BEM.aT.b2','BEM.aT.b3','BEM.aT_new.b1','BEM.aT_new.b2','BEM.aT_new.b3',...
                  'BEM.V_Res.b1','BEM.V_Res.b2','BEM.V_Res.b3','BEM.phi.b1','BEM.phi.b2','BEM.phi.b3',...
                  'BEM.alpha.b1','BEM.alpha.b2','BEM.alpha.b3','BEM.CL.b1','BEM.CL.b2','BEM.CL.b3',...
                  'BEM.CD.b1','BEM.CD.b2','BEM.CD.b3','BEM.sig',...
                  'BEM.CN.b1','BEM.CN.b2','BEM.CN.b3','BEM.CT.b1','BEM.CT.b2','BEM.CT.b3',...
                  'BEM.f.b1','BEM.f.b2','BEM.f.b3','BEM.F.b1','BEM.F.b2','BEM.F.b3',...
                  'BEM.iter_count.b1','BEM.iter_count.b2','BEM.iter_count.b3','BEM.FlapdF.b1','BEM.FlapdF.b2','BEM.FlapdF.b3',...
                  'BEM.Rad_dF.b1','BEM.Rad_dF.b2','BEM.Rad_dF.b3','BEM.EdgedF.b1','BEM.EdgedF.b2','BEM.EdgedF.b3',...
                  'BEM.FlapdM.b1','BEM.FlapdM.b2','BEM.FlapdM.b3','BEM.EdgedM.b1','BEM.EdgedM.b2','BEM.EdgedM.b3',...
                  'BEM.EdgedMC.b1','BEM.EdgedMC.b2','BEM.EdgedMC.b3','BEM.dT','BEM.Ct','BEM.Q_L.b1','BEM.Q_L.b2','BEM.Q_L.b3',...
                  'BEM.FlapdF_max','BEM.EdgedF_max','BEM.FlapdM_max','BEM.EdgedM_max','BEM.EdgedMC_all'); 
end