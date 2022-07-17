 function [output_details, output, BEM] = core_bem(General, op_pts, BLD)
    % Core steady state BEM calculations
    % Comments
    % 
    % 
    
    % Operating data organising
    wsp = op_pts.wsp;
    shear = zeros(length(wsp),1);
    pitch = op_pts.pitch;
    rpm = op_pts.rpm;
    
    if all(shear >= 0), bld_calcs = 1; else bld_calcs =3; end
    
    BEM = initialize_structs(); % Pre-initializing structure for MATLAB Codegen
              
    BEM.r = BLD.r; % Blade stations
    BEM.C = BLD.C; % Chord distribution
    BEM.t_C = BLD.t_C; % t/C distribution
    BEM.AeroTwist = BLD.AeroTwist; % Aerodynamic (mold) Twist distribution [deg]
    BEM.preflap = BLD.preflap; % Pre-bend distribution (m)
    BEM.PreBendTwist = calc_twist(BEM.r, BEM.preflap); % Pre-bend angle w.r.t root [deg]
    
    BEM.Omega_r = (repmat(rpm',length(BEM.r),1)*2*pi/60).*repmat(BEM.r,1,length(rpm)); % Velocity seen by each of the section due to rotation [m/s]
    BEM.sig = BEM.C*General.N./(2*pi.*BEM.r); % Solidity ratio
    
    b = fieldnames(BEM.aA); % Getting the blade struct field names
    for id=1:numel(b) % Loop for each of the blades
        if id <= bld_calcs % Condition for whether to skip 2nd and 3rd blade calculation
            BEM.aA.(b{id}) = zeros(length(BEM.r),size(wsp,1)); % Initializing axial induction factors 
            BEM.aT.(b{id}) = zeros(length(BEM.r),size(wsp,1)); % Initializing tangential induction factors 
            BEM.V_Res.(b{id}) = sqrt(repmat((wsp.^2)',length(BEM.r),1)+BEM.Omega_r.^2); % Resultant velocity seen by each of the section [m/s]
            BEM.phi.(b{id}) = atan2d((repmat(wsp',length(BEM.r),1).*(1-BEM.aA.(b{id}))),(BEM.Omega_r.*(1+BEM.aT.(b{id})))); % Calculating flow angle PHI [deg]
            BEM.alpha.(b{id}) = BEM.phi.(b{id})-repmat(BEM.AeroTwist,1,length(wsp))-repmat(pitch',length(BEM.r),1); % Angle of attack seen by each section [deg]
            [BEM.CL.(b{id}), BEM.CD.(b{id})] = CL_CD_vs_alpha(BEM.t_C, BLD.pro_t_C, BLD.pro_AoA, BLD.pro_cL, BLD.pro_cD, BEM.alpha.(b{id})); % CL, CD coefficients
            BEM.CN.(b{id}) = BEM.CL.(b{id}).*cosd(BEM.phi.(b{id}))+BEM.CD.(b{id}).*sind(BEM.phi.(b{id})); % Sectional Normal force coefficient
            BEM.CT.(b{id}) = BEM.CL.(b{id}).*sind(BEM.phi.(b{id}))-BEM.CD.(b{id}).*cosd(BEM.phi.(b{id})); % Sectional Tangential force coefficient

            if General.induction % Condition for doing induction calculations
                BEM.aA_new.(b{id}) = zeros(length(BEM.r),size(wsp,1)); % Initializing axial induction factors 
                BEM.aT_new.(b{id}) = zeros(length(BEM.r),size(wsp,1)); % Initializing tangential induction factors 
                BEM.f.(b{id}) = zeros(length(BEM.r),size(wsp,1)); % Initializing f matrix
                BEM.F.(b{id}) = zeros(length(BEM.r),size(wsp,1)); % Initializing f matrix
                BEM.iter_count.(b{id}) = zeros(length(BEM.r),size(wsp,1)); % Initializing f matrix
                for i=1:length(BEM.r) % Loop for different sections of the blade
                    for j=1:size(wsp,1) % Loop for different operating set points
                        BEM.aA_new.(b{id})(i,j) = 1./((4.*sind(BEM.phi.(b{id})(i,j)).*sind(BEM.phi.(b{id})(i,j))./(BEM.sig(i).*BEM.CN.(b{id})(i,j)))+1); % Calculating new axial induction factors
                        BEM.aT_new.(b{id})(i,j) = 1./((4.*sind(BEM.phi.(b{id})(i,j)).*cosd(BEM.phi.(b{id})(i,j))./(BEM.sig(i).*BEM.CT.(b{id})(i,j)))-1); % Calculating new axial induction factors
                        iter=1;
                        while (abs(BEM.aA_new.(b{id})(i,j) - BEM.aA.(b{id})(i,j)) > 0.0001 || abs(BEM.aT_new.(b{id})(i,j) - BEM.aT.(b{id})(i,j)) > 0.0001) && iter < 1000 % Induction iteration
                            if BEM.aA_new.(b{id})(i,j) > 1.5, BEM.aA.(b{id})(i,j) = 1.5; elseif BEM.aA_new.(b{id})(i,j) < -1, BEM.aA.(b{id})(i,j) = -1; else, BEM.aA.(b{id})(i,j) = BEM.aA_new.(b{id})(i,j);end
                            if BEM.aT_new.(b{id})(i,j) > 1,   BEM.aT.(b{id})(i,j) = 1;   elseif BEM.aT_new.(b{id})(i,j) < -1, BEM.aT.(b{id})(i,j) = -1; else, BEM.aT.(b{id})(i,j) = BEM.aT_new.(b{id})(i,j);end
                            BEM.phi.(b{id})(i,j) = atan2d((wsp(j)'.*(1-BEM.aA.(b{id})(i,j))),(BEM.Omega_r(i,j).*(1+BEM.aT.(b{id})(i,j)))); % Calculating flow angle PHI [deg]
                            if General.tip_loss
                                BEM.f.(b{id})(i,j) = (General.N/2).*(BEM.r(end)-BEM.r(i))./(BEM.r(i).*sind(BEM.phi.(b{id})(i,j)));
                                if BEM.f.(b{id})(i,j)<0 % Reverse flow situation, Aero twist needs to be adjusted
                                    BEM.f.(b{id})(i,j) = 999;
                                end
                                BEM.F.(b{id})(i,j)=(2/pi).*acos(exp(-BEM.f.(b{id})(i,j)));
                            else
                                BEM.F.(b{id})(i,j)=1;
                            end
                            BEM.alpha.(b{id})(i,j) = BEM.phi.(b{id})(i,j)-BEM.AeroTwist(i)-pitch(j)'; % Angle of attack seen by each section [deg]
                            if BEM.alpha.(b{id})(i,j) < -180 || BEM.alpha.(b{id})(i,j) > 180
                                BEM.alpha.(b{id})(i,j) = mod(BEM.alpha.(b{id})(i,j), 180); % Check here once
                            end
                            [BEM.CL.(b{id})(i,j), BEM.CD.(b{id})(i,j)] = CL_CD_vs_alpha(BEM.t_C(i), BLD.pro_t_C, BLD.pro_AoA, BLD.pro_cL, BLD.pro_cD, BEM.alpha.(b{id})(i,j)); % CL, CD coefficients
                            BEM.CN.(b{id})(i,j) = BEM.CL.(b{id})(i,j).*cosd(BEM.phi.(b{id})(i,j))+BEM.CD.(b{id})(i,j).*sind(BEM.phi.(b{id})(i,j)); % Sectional Normal force coefficient (Thrust direction force coefficient)
                            BEM.CT.(b{id})(i,j) = BEM.CL.(b{id})(i,j).*sind(BEM.phi.(b{id})(i,j))-BEM.CD.(b{id})(i,j).*cosd(BEM.phi.(b{id})(i,j)); % Sectional Tangential force coefficient
                            if General.highCT
                                method = General.highCT; % Method = 1 or 2
                                BEM.aA_new.(b{id})(i,j) = real(calc_ind_using_high_CT_approx(method, BEM.aA.(b{id})(i,j), BEM.F.(b{id})(i,j), BEM.phi.(b{id})(i,j), BEM.sig(i), BEM.CN.(b{id})(i,j)));
                            else
                                BEM.aA_new.(b{id})(i,j) = 1./((4.*BEM.F.(b{id})(i,j).*sind(BEM.phi.(b{id})(i,j)).*sind(BEM.phi.(b{id})(i,j))./(BEM.sig(i).*BEM.CN.(b{id})(i,j)))+1); % Calculating new axial induction factors
                            end
                            BEM.aT_new.(b{id})(i,j) = 1./((4.*BEM.F.(b{id})(i,j).*sind(BEM.phi.(b{id})(i,j)).*cosd(BEM.phi.(b{id})(i,j))./(BEM.sig(i).*BEM.CT.(b{id})(i,j)))-1); % Calculating new axial induction factors
                            RF=0.50; % Relaxation Factor as per https://onlinelibrary.wiley.com/doi/full/10.1002/ese3.945
                            BEM.aA_new.(b{id})(i,j) = RF*BEM.aA_new.(b{id})(i,j)+(1-RF)*BEM.aA.(b{id})(i,j);
                            BEM.aT_new.(b{id})(i,j) = RF*BEM.aT_new.(b{id})(i,j)+(1-RF)*BEM.aT.(b{id})(i,j);
                            iter=iter+1;
                        end
                        BEM.iter_count.(b{id})(i,j)=iter; % Induction iteration counter
                    end
                end
            end
        end
    end
        
    
    %% Sectional Calculations
    A = pi.*BEM.r(end).^2; % Actual Rotor area
    Q_G = 0.5.*General.rho.*wsp'.*wsp'; % Dynamic pressure
    
    % Sectional loads in Blade root coords
    for id=1:numel(b)
        if id <= bld_calcs
            BEM.Q_L.(b{id}) = 0.5*General.rho.*BEM.V_Res.(b{id}).*BEM.V_Res.(b{id}); % Local Dynamic pressure for each station
            BEM.FlapdF.(b{id}) = BEM.Q_L.(b{id}).*repmat(BEM.C,1,length(wsp)).*BEM.CN.(b{id}).*repmat(cosd(BEM.PreBendTwist),1,length(wsp)); % Local Thrust/Normal force [N], Check this
            BEM.Rad_dF.(b{id}) = BEM.Q_L.(b{id}).*repmat(BEM.C,1,length(wsp)).*BEM.CN.(b{id}).*repmat(sind(BEM.PreBendTwist),1,length(wsp)); % Local Radial force [N], Check this
            BEM.EdgedF.(b{id}) = BEM.Q_L.(b{id}).*repmat(BEM.C,1,length(wsp)).*BEM.CT.(b{id}); % Local Edge force [N]
            BEM.FlapdM.(b{id}) = BEM.FlapdF.(b{id}).*repmat(BEM.r-BEM.r(1),1,length(wsp))+BEM.Rad_dF.(b{id}).*repmat((-BEM.preflap),1,length(wsp)); % Local Flap Moment [Nm] at Hub end (or Blade root)
            BEM.EdgedM.(b{id}) = BEM.EdgedF.(b{id}).*repmat(BEM.r-BEM.r(1),1,length(wsp)); % Local Edge Moment [Nm] at Hub end (or Blade root)
        end
    end
    
    if bld_calcs == 1
        BEM.FlapdF_max = BEM.FlapdF.b1;
        BEM.EdgedF_max = BEM.EdgedF.b1;
        BEM.FlapdM_max = BEM.FlapdM.b1;
        BEM.EdgedM_max = BEM.EdgedM.b1;
        BEM.EdgedMC_all = General.N*BEM.EdgedF.b1.*repmat(BEM.r,1,length(wsp));
        BEM.dT = General.N*BEM.FlapdF.b1;
    else
        BEM.FlapdF_max = max(max(BEM.FlapdF.b1, BEM.FlapdF.b2), BEM.FlapdF.b3);
        BEM.EdgedF_max = max(max(BEM.EdgedF.b1, BEM.EdgedF.b2), BEM.EdgedF.b3);
        BEM.FlapdM_max = max(max(BEM.FlapdM.b1, BEM.FlapdM.b2), BEM.FlapdM.b3);
        BEM.EdgedM_max = max(max(BEM.EdgedM.b1, BEM.EdgedM.b2), BEM.EdgedM.b3);
        BEM.EdgedMC_all = (BEM.EdgedF.b1 + BEM.EdgedF.b2 + BEM.EdgedF.b3).*repmat(BEM.r,1,length(wsp)); % Local Edge Moment [Nm] (w.r.t rotor center) in blade root coords
        BEM.dT = BEM.FlapdF.b1+BEM.FlapdF.b2+BEM.FlapdF.b3; % Full Rotor loads
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