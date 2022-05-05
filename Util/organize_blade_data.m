function Blade = organize_blade_data(sec_r, sec_C, sec_t_C, sec_twist, pc)
    % This function organizes all the blade structural and aerodynamic data
    % in a structure for better data handling
    
    % sec_r [m] is is the blade section table (m x 1 vector)
    % sec_C [m] is the Chord corresponding to each of the sec_r (m x 1 vector)
    % sec_t_C is the thickness to chord ratio corresponding to each of the sec_r (m x 1 vector)
    % sec_twist [deg] is is the Aerodynamic twist at each of the sec_r (m x 1 vector)
    % preflap [m] is the section wise offset of each of the sections away from tower (m x 1 vector)
    % profile_t_C is the K x 1 thickneses available
    % profile_AoA is the L x 1 angles of attack
    % profile_cL is L x K cL values corresponding to each angles of attack and t_C
    % profile_cD is L x K cD values corresponding to each angles of attack and t_C
    
    % Input data organising, Making sure values are in column vector
    sec_r = sec_r(:);
    sec_C = sec_C(:);
    sec_t_C = sec_t_C(:);
    sec_twist = sec_twist(:);
    
    ae = readmatrix('G:\BEM\BEM\Data\NREL5MWRefTurb_v50\data\NREL_5MW_ae.txt');
    
    % Organizing the blade data in a structure
    Blade.sec_r = sec_r;
    Blade.sec_C = sec_C;
    Blade.sec_t_C = sec_t_C;
    Blade.sec_twist = sec_twist;
    Blade.preflap = zeros(length(sec_r),1);
    
    [Blade.pro_t_C, Blade.pro_AoA, Blade.pro_cL, Blade.pro_cD] = read_pc_file(pc);
end