%% MAIN PERCOLATION-DIFFUSION
% Computes the diffusional re-equilibration of m trace elements in a solid matrix
% of n spherical minerals percolated by a fluid in a 1D column
% Flexible MATLAB implementation of a 1D percolation-diffusion model computing the diffusional re-equilibration of m trace elements in a solid matrix of n spherical phases percolated by a liquid

% Modified from the original TPMCRT code developed by Beñat Oliveira Bravo
% By Romain Tilhac, May 2024

    % Development notes:
    % When coupled
    % check initializeDynamics for non-uniform input arrays
    % check why density are set to 1 in rho_nTP
    % check why modal proportion are passed to vol. fraction
    % check 2D comments from original version?

%% LOAD INPUTS

    cd(fileparts(mfilename('fullpath'))); % ensure that the filepath is open
    restoredefaultpath
    addpath(fullfile(pwd));% add model path
    addpath(fullfile(pwd,'utils'));% add utils pathf
    
    clc; % clear command
    close all; % close opened figures
    fclose('all'); % close opened files
    clear; % clear workspace
    closevar; %close opened variables

    % ====================load M3G======================
    if ismac
        %matlab_path='/Users/romaintilhac/Documents/Research/computing/MATLAB_local/';
        matlab_path='/Users/romaintilhac/Nextcloud/MATLAB/';
    elseif isunix
        matlab_path='/home/romaintilhac/Nextcloud/MATLAB/';
    end
    addpath(fullfile(matlab_path, 'M3G/pre-release_1.0'));
    model_type = "percolation-diffusion";
    setup_name="default";
    M3G = initializeM3G(model_type, fullfile(pwd,'setups',setup_name));
    disp('M3G package succesfully loaded.');
    % =================================================

    % unpack M3G    
    lists=M3G.lists;
    settings=M3G.settings;
    options=M3G.options;
    parameters=M3G.parameters;
    compositions=M3G.compositions;
    coefficients=M3G.coefficients;
    Kd=M3G.Kd;
    diffusivities=M3G.diffusivities;

    % read lists
    TP_list = lists.TP; nTP = lists.nTP;
    TE_list = lists.TE; nTE = lists.nTE;

    % read options
    fix_radi = options.fix_radi;
    Eu_diff = options.Eu_diff;
    syntheticCpx = options.syntheticCpx;

    % read settings
    elemTypeT = settings.elemTypeT;
    nent = settings.nent;
    nxt = settings.nxt;
    nyt = settings.nyt;
    mxelem = settings.mxelem;
    myelem = settings.myelem;
    markmove = settings.markmove;
    nc = settings.nc;
    exp_factor = settings.exp_factor;
    n_Tp0 = settings.n_Tp0;
    tol_NTE = settings.tol_NTE;
    tol_TE = settings.tol_TE;
    tol_TE_abs = settings.tol_TE_abs;
    max_iter_TE = settings.max_iter_TE;
    iter_negative_max = settings.iter_negative_max;
    damping0 = settings.damping0;
    damping1 = settings.damping1;
    save_interval = settings.save_interval;

    % read scalar parameters
    time_end = parameters.time_end;
    x1 = parameters.x1;
    x2 = parameters.x2;
    y1 = parameters.y1;
    y2 = parameters.y2;
    fix_poros = parameters.fix_poros;
    fix_P = parameters.fix_P;
    fix_T = parameters.fix_T;
    fix_Vy = parameters.fix_Vy;
    r_Eu = parameters.r_Eu;

    % initialization
    nTP_vol = zeros(1, nTP);
    diff_nTE = ones(1,nTE);
    ini_radi=zeros(1, nTP); 
    Kd_coeff=zeros(nTP, nTE);
    D_coeff=zeros(nTP, nTE);
    E_coeff=zeros(nTP, nTE);
    V_coeff=zeros(nTP, nTE);
    
    % read compositions & nTp arrays
    TE_input_solid = table2array(compositions.TE.Solid);
    TE_input_liquid = table2array(compositions.TE.Liquid);
    for iTP=1:nTP
        TP_str=TP_list(iTP);
        diffT_nTP(iTP)=options.("diffT_TE_Solid_"+TP_str);
        diffP_nTP(iTP)=options.("diffP_TE_Solid_"+TP_str)*diffT_nTP(iTP);
        ini_radi(iTP)=parameters.("ini_radi_Solid_"+TP_str);
        nTP_vol(iTP)=coefficients.mode_Solid.(TP_str);
        Kd_coeff(iTP, :)=table2array(Kd.mineralKd.TE.Solid.(TP_str));
        D_coeff(iTP, :)=table2array(diffusivities.TE.D.Solid.("D_"+TP_str));
        if diffT_nTP(iTP)==1; E_coeff(iTP, :)=table2array(diffusivities.TE.E.Solid.("E_"+TP_str))*diffT_nTP(iTP); end
        if diffP_nTP(iTP)==1; V_coeff(iTP, :)=table2array(diffusivities.TE.V.Solid.("V_"+TP_str))*diffP_nTP(iTP); end
    end
    nTP_vol(isnan(nTP_vol)) = 0;
    diff_nTP = double(nTP_vol>0);
    diffT_nTP = double(nTP_vol>0);
    diffP_nTP = double(nTP_vol>0);

% % read input file & set default values
% 
% filename='input.xlsx';
% 
% TE_input_solid=readtable(filename).solid';
% TE_input_liquid=readtable(filename).liquid';
% 
% Kd_coeff=zeros(nTP, nTE);
% D_coeff=zeros(nTP, nTE);
% E_coeff=zeros(nTP, nTE);
% V_coeff=zeros(nTP, nTE);
% diff_nTE = ones(1,nTE);
% 
% for iTP=1:nTP
%     mineral=TP_list(iTP);
%     Kd.(mineral)=readtable(filename).(strcat('Kd_',mineral))';
%     D.(mineral)=readtable(filename).(strcat('D_',mineral))';
%     E.(mineral)=readtable(filename).(strcat('E_',mineral))';
%     V.(mineral)=readtable(filename).(strcat('V_',mineral))';
%     Kd_coeff(iTP,:) = Kd.(mineral);
%     D_coeff(iTP,:) = D.(mineral);
%     E_coeff(iTP,:) = E.(mineral)*diffT_nTP(iTP);
%     V_coeff(iTP,:) = V.(mineral)*diffP_nTP(iTP);
% end
% 
% % TP and TE must be consistent with input file
% TP_list= split("Cpx Cpx"); nTP = length(TP_list);
% TE_list=split("La Ce Pr Nd Sm Eu Gd Tb Dy Ho Er Tm Yb Lu"); nTE = length(TE_list);
% 
% nTP_vol =                                     [0.9    0.1];  % modal proportion [vol]
% diff_nTP  =  double(nTP_vol>0);
% ini_radi =                                      [0.01   1];  % grain size
% diffT_nTP =  diff_nTP.*               [1      1];  % T dependency of the diffusivity
% diffP_nTP = diffT_nTP.*              [1      1];  % P dependency of the diffusivity
% 
% fix_radi            = 0;     % 1, fix grain size to initial sizeEu_diff
% Eu_diff = 1; % to use Eu diffusivities in Cpx based on Sr diffusivities from Sneeringer et al. 1984 (only used if Eu and Cpx are present)
% syntheticCpx = 0; % to use 1 synthetic diopside (conservative), 0 to use natural diopside (only used if Eu_diff == 1)
% 
% % scalar parameters
% time_end            = 0.005; % duration [My]
% x1                  = 0;     % domain X first vertice [m]
% x2                  = 300;   % domain X second vertice [m]     
% y1                  = 0;     % domain Y first vertice [m]
% y2                  = 1000;  % domain Y second vertice  [m]
% fix_poros           = 0.01;  % constant melt porosity (vol. fraction)
% fix_P               = 1.2;   % constant pressure [GPa]
% fix_T               = 1200;  % constant temperature [ºC]
% fix_Vy              = 5;     % constant melt velocity [cm/y]
% r_Eu = 0.1; %Eu2+/Eu3+ proportion (only used if options.Eu_diff == 1)
% 
% % computational settings
% elemTypeT           = 1;     % Quads (Eulerian mesh for temperature)
% nent                = 4;     % Number of nodes per element (velocity)
% nxt                 = 10;    % Number of nodes -1 in x direction
% nyt                 = 10;    % Number of nodes -1 in y direction
% mxelem              = 5;     % Number of Lagrangian markers in x direction per element WARNING unused
% myelem              = 5;     % Number of Lagrangian markers in y direction per element
% markmove            = 4;     % Moving Lagrangian markers: 0 = not moving at all, 1 = simple 1-st order advection, 4 = 4-th order in space Runge-Kutta
% nc	                = 51;	 % number of nodes for the diffusion profiles
% exp_factor          = 0.3;   % exponential factor for non-uniform mesh generation
% n_Tp0               = 5E5;   % default number of particles per unit volume
% tol_NTE	            = 1E-8;	 % numerical tolerance for non-trivial existence of phase (mass or vol. fraction)
% tol_TE              = 1E-9;  % numerical tolerance for the relative error in trace element concentrations (relative)
% tol_TE_abs          = 5E-9;  % numerical tolerance for convergence of trace element concentrations (absolute mass fraction)
% max_iter_TE         = 40;    % maximum number of iterations
% iter_negative_max   = 30;    % number of correction iterations to fix negative fluid compositions 
% damping0            = 0.5;   % initial damping factor for TE iteration
% damping1            = 1;     % damping factor for TE iteration
% save_interval       = 1;     % saving interval: 1, every timestep

% unit conversions
time_end = time_end*1e6*365*24*60*60;       % [My]   to [s]
fix_Vy = fix_Vy/100/365/24/60/60;           % [cm/y] to [m/s]
ini_radi = ini_radi.*1E-3;                  % [mm]   to [m]

melt_idx = nTP + 1;
   
%% SETUP PROBLEM

    % PROBLEM STATEMENT 
        xsize = x2 - x1;
        ysize = y2 - y1;

    % EULERIAN MESH
        % Temperature
        nnt = (nxt+1)*(nyt+1);       % Number of total temperature nodes
        xstpt = xsize/nxt;           % Horizontal grid step
        ystpt = ysize/nyt;           % Vertical grid step
            
        % Create the meshes for temperature
        [XT,TT] = createVelocityMesh(elemTypeT,nent,x1,x2,y1,y2,nxt,nyt);

        % Making vectors for nodal points positions (basic nodes)
            gridyt   =   (0:ystpt:ysize);   % Vertical
            gridxt   =   (0:xstpt:xsize);   % Horizontal

    % LAGRANGIAN PARTICLES
        % Create the particles
        [ME_solid, ME_fluid]  = createParticles_Elements(mxelem,myelem,nxt,nyt,xsize,ysize,nTE,nTP);

    [NR_solid,NR_fluid,NGamma,NRho0,NRho,NRho_tp,NRho_tp0,NV_mat0,NV_mat,NP_mat,NDiv0,NDiv,NPhi0,NPhi,NTp_wt0,NTp_wt,NTp_vol0,NTp_vol,NT] = ...
    initializeDynamics('default.mat',nxt,nyt);

    timestep = ystpt/fix_Vy/myelem;
    ntimesteps = round(time_end/timestep);
    save_list = 1:save_interval:ntimesteps;

    rho_nTP = diff_nTP.*ones(1, nTP);
    rho_nPH   = [nTP_vol*rho_nTP' 1];
    Rho_ave   = [(1-fix_poros) fix_poros]*rho_nPH';
    input_TP_v = [nTP_vol/sum(nTP_vol)*(1-fix_poros) fix_poros];
    input_TP_w = input_TP_v.*[rho_nTP rho_nPH(2)]./Rho_ave;
    
    n_Tp = repmat(n_Tp0, 1, nTP);
    n_Tp(logical(diff_nTP)) = (nTP_vol(logical(diff_nTP))/sum(nTP_vol)*(1-fix_poros))./((4/3)*pi.*ini_radi(logical(diff_nTP)).^3);
    ME_solid(:,3+nTE+1:end) = repmat(n_Tp,size(ME_solid,1),1);
    
    fix_radi  = fix_radi*diff_nTP.*ini_radi; 
    fix_poros = [1-fix_poros fix_poros];    
    fix_press = [fix_P fix_P];
    fix_velo  = [0 0 0 fix_Vy];

    NR_solid = 0*NR_solid;
    NR_fluid = 0*NR_fluid;
    NGamma = 0*NGamma;
    NRho0 = ones(size(NRho0));
    NRho =  ones(size(NRho));
    NRho_tp = repmat(rho_nTP,size(NRho_tp,1),1);
    NRho_tp0 = repmat(rho_nTP,size(NRho_tp0,1),1);
    NV_mat0 = repmat(fix_velo,size(NV_mat0,1),1);
    NV_mat = repmat(fix_velo,size(NV_mat,1),1);
    NP_mat = repmat(fix_press,size(NP_mat,1),1);
    NDiv0 = 0*NDiv0;
    NDiv = 0*NDiv;
    NPhi0 = repmat(fix_poros,size(NPhi0,1),1);
    NPhi = repmat(fix_poros,size(NPhi,1),1);
    NTp_wt0 = repmat(input_TP_w,size(NTp_wt0,1),1);
    NTp_wt = repmat(input_TP_w,size(NTp_wt,1),1);
    NTp_vol0 = repmat(input_TP_v,size(NTp_vol0,1),1);
    NTp_vol = repmat(input_TP_v,size(NTp_vol,1),1);
    NT = repmat(fix_T,size(NT,1),1);

    fprintf('Initialization completed - %d timesteps \n', ntimesteps)
    keyboard

%% RUN MODEL
        
    iter = 1;  save_iter = 0;
    time = 0; timeMy_mat = [];
    
    while time<time_end
         
        if iter==1
    
            %% FIRST TIMESTEP
    
            % Flexible per-phase cell Storage
    
            ME_tp      = cell(1, nTP);    % Trace elements per phase
            ME_eq_tp   = cell(1, nTP);
            Kd_tp      = cell(1, nTP);    % Partition coefficients per phase
            D_tp       = cell(1, nTP);    % Diffusion coefficients per phase
            E_tp       = cell(1, nTP);    % Activation energies per phase
            V_tp       = cell(1, nTP);    % Activation volumes per phase
            Mdiff_tp   = cell(1, nTP);
            Gamma_tp   = cell(1, nTP);    % Gamma for each phase
            MBC_tp     = cell(1, nTP);    % Boundary conditions per phase
            MBC0_tp    = cell(1, nTP);
            Mat_aux_tp = cell(1, nTP);    % For TE profiles
            Rradi_0    = cell(1, nTP);
            TEprofile_tp= cell(nTP, nTE);
            Error_tp = cell(1, nTP);
            MBC_iter_tp = cell(1, nTP);
            ME_noChange_tp = cell(1, nTP);
            TE_noChange_tp = cell(1, nTP);
            nSolid = size(ME_solid, 1);
            nFluid = size(ME_fluid, 1);
            auxE_tp = cell(1, nTP+1); % One extra for melt
            ME_solid_RhoPhi_tp = cell(1, nTP+1);
            ME_fluid_RhoPhi_tp = cell(1, nTP+1);
            ME_solid_w_tp = cell(1, nTP+1);
    
            for iTP = 1:nTP
                ME_tp{iTP}    = zeros(nSolid, nTE);
                ME_eq_tp{iTP} = zeros(nSolid, nTE);
                Kd_tp{iTP}    = repmat(Kd_coeff(iTP,:), nSolid, 1); % Rows # particles, columns TE
                D_tp{iTP}     = repmat(D_coeff(iTP,:), nSolid, 1);
                E_tp{iTP}     = repmat(E_coeff(iTP,:), nSolid, 1);
                V_tp{iTP}     = repmat(V_coeff(iTP,:), nSolid, 1);
                Mdiff_tp{iTP} = zeros(nSolid, nTE);
                Gamma_tp{iTP} = zeros(nSolid, nTE);
                MBC_tp{iTP}   = zeros(nSolid, nTE);
                MBC0_tp{iTP}  = zeros(nSolid, nTE);
            end
        
            % Flexible matrices (TP or nTP+1 cols)
    
            ME_fluid_TE = zeros(nFluid, nTE);
            ME_solid_Rho_TP  = zeros(nSolid, nTP+1);
            ME_fluid_Rho_TP  = zeros(nFluid, nTP+1);
            
            for iTP = 1:nTP
                ME_solid_Rho_TP(:, iTP) = interp1(XT(1:nxt+1:end,2), NRho_tp0(1:nxt+1:end,iTP), ME_solid(:,2));
                ME_fluid_Rho_TP(:, iTP) = interp1(XT(1:nxt+1:end,2), NRho_tp0(1:nxt+1:end,iTP), ME_fluid(:,2));
            end
            NRho0(NRho0(:,2)==0,2) = NRho0(NRho0(:,2)==0,1);
            ME_solid_Rho_TP(:, melt_idx) = interp1(XT(1:nxt+1:end,2), NRho0(1:nxt+1:end,2), ME_solid(:,2));
            ME_fluid_Rho_TP(:, melt_idx) = interp1(XT(1:nxt+1:end,2), NRho0(1:nxt+1:end,2), ME_fluid(:,2));
        
            ME_solid_Rho_TP0 = ME_solid_Rho_TP;
    
            ME_solid_TP_w0 = zeros(nSolid, nTP+1);
            ME_fluid_TP_w0 = zeros(nFluid, nTP+1);
            ME_solid_TP_v0 = zeros(nSolid, nTP+1);
            ME_fluid_TP_v0 = zeros(nFluid, nTP+1);
            for iTP = 1:nTP+1
                ME_solid_TP_w0(:, iTP) = interp1(XT(1:nxt+1:end,2), NTp_wt0(1:nxt+1:end,iTP), ME_solid(:,2));
                ME_fluid_TP_w0(:, iTP) = interp1(XT(1:nxt+1:end,2), NTp_wt0(1:nxt+1:end,iTP), ME_fluid(:,2));
                ME_solid_TP_v0(:, iTP) = interp1(XT(1:nxt+1:end,2), NTp_vol0(1:nxt+1:end,iTP), ME_solid(:,2));
                ME_fluid_TP_v0(:, iTP) = interp1(XT(1:nxt+1:end,2), NTp_vol0(1:nxt+1:end,iTP), ME_fluid(:,2));
            end
            ME_solid_TP_w0 = ME_solid_TP_w0 ./ repmat(sum(ME_solid_TP_w0,2),1, size(ME_solid_TP_w0,2));
            ME_fluid_TP_w0 = ME_fluid_TP_w0 ./ repmat(sum(ME_fluid_TP_w0,2),1, size(ME_fluid_TP_w0,2));
            ME_solid_TP_v0 = ME_solid_TP_v0 ./ repmat(sum(ME_solid_TP_v0,2),1, size(ME_solid_TP_v0,2));
            ME_fluid_TP_v0 = ME_fluid_TP_v0 ./ repmat(sum(ME_fluid_TP_v0,2),1, size(ME_fluid_TP_v0,2));
            
            for iTP = 1:nTP
                Kd_tp{iTP} = repmat(Kd_coeff(iTP,:), nSolid, 1);
            end
    
            % Initial TE (ppm) concentration in TP and Melt 
            ME_TE = repmat(TE_input_solid, nSolid, 1);
    
            ME_solid_w0 = zeros(nSolid, nTE, nTP);
            for iTP = 1:nTP+1
                ME_solid_w0(:,:,iTP) = repmat(ME_solid_TP_w0(:, iTP), 1, nTE);
            end
            ME_solid_w0_solid = sum(ME_solid_w0(:,:,1:nTP), 3); % normalize to 1 every TP contribution
            ME_solid_w0_total = sum(ME_solid_w0, 3);% normalize to 1 every TP contribution
        
            ME_solid(:,4:3+nTE) = ME_TE; % input is the solid
            ME_solid0 = ME_solid;
    
            auxE = zeros(nSolid, nTE, nTP);
            denom = zeros(nSolid, nTE);
            for iTP = 1:nTP
                auxE(:,:,iTP) = repmat(ME_solid_TP_v0(:, iTP) > 0, 1, nTE);
                denom = denom + auxE(:,:,iTP) .* Kd_tp{iTP} .* ME_solid_w0(:,:,iTP) ./ ME_solid_w0_solid;
            end
            ME_melt_aux = ME_solid(:,4:3+nTE) ./ denom; % [melt] = [solid]/Kd_bulk    (normalized to solid)
            ME_melt = zeros(nSolid, nTE);
            meltMask = ME_solid_TP_v0(:, melt_idx) > 0;
            ME_melt(meltMask, :) = ME_melt_aux(meltMask, :); % TE mass / fluid mass 
            ME_melt_eq = ME_melt_aux;
        
            for TE_pos = 1:nTE
                ME_fluid(:,3+TE_pos) = interp1(ME_solid(:,2), ME_melt(:,TE_pos), ME_fluid(:,2), 'nearest', 'extrap');
            end
    
            % TE mass / TP mass 
            for iTP = 1:nTP
                ME_eq_tp{iTP} = ME_melt_aux .* auxE(:,:,iTP) .* Kd_tp{iTP};
                phaseMask = ME_solid_TP_v0(:,iTP) > 0; % only where phase is present
                ME_tp{iTP}(phaseMask,:) = ME_eq_tp{iTP}(phaseMask,:);
            end
        
            % Prepare BC for next iteration
            for iTP = 1:nTP
                MBC0_tp{iTP} = ME_tp{iTP};
                MBC0_tp{iTP}(meltMask,:) = ME_melt(meltMask,:) .* Kd_tp{iTP}(meltMask,:);
                MBC_tp{iTP} = MBC0_tp{iTP};
            end
            
            % TE profile arrays
            for iTP = 1:nTP
                for TE_pos = 1:nTE
                    TEprofile_tp{iTP,TE_pos} = repmat(ME_tp{iTP}(:,TE_pos)', nc, 1)';
                end
                Mat_aux_tp{iTP} = cell2mat(TEprofile_tp(iTP,:));
                Mat_order = reshape(Mat_aux_tp{iTP}', nc, [])';
                ME_tp{iTP} = mat2cell(Mat_order, nTE*ones(1, nSolid), nc);
            end
            for iTP = 1:nTP
                Mat_order = reshape(Mat_aux_tp{iTP}', nc, [])';
                ME_tp{iTP} = mat2cell(Mat_order, nTE*ones(1, nSolid), nc);
            end
            ME0_tp = ME_tp;
    
            % Diffusion
    
            % P-T in particles
            ME_T = interp1(XT(1:nxt+1:end,2), NT(1:nxt+1:end,1), ME_solid(:,2)) + constants.TK;
            ME_P = interp1(XT(1:nxt+1:end,2), NP_mat(1:nxt+1:end,1), ME_solid(:,2)) * constants.PGPa;
            
            % Partition coefficients & diffusivites
            for iTP = 1:nTP
                Dmat = D_tp{iTP}; 
                Emat = E_tp{iTP};
                Vmat = V_tp{iTP};
                Mdiff_tp{iTP} = double(Dmat .* exp((-Emat + repmat(ME_P,1,nTE).*Vmat) ./ (constants.R.*repmat(ME_T,1,nTE))));
            end
            
            % Calculate Eu diffusivity in Cpx
            if Eu_diff
                iCpx =  find(TP_list=='Cpx');
                if isempty(iCpx); error('Cpx not found in TP_list.'); end
                if isscalar(iCpx)
                    Mdiff_tp{iCpx}  = calcEuCpxdiffusivity(Mdiff_tp{iCpx}, ME_T, ME_P, TE_list, r_Eu, syntheticCpx);
                else
                    for i = 1:length(iCpx)
                        Mdiff_tp{iCpx(i)}  = calcEuCpxdiffusivity(Mdiff_tp{iCpx(i)}, ME_T, ME_P, TE_list, r_Eu, syntheticCpx);
                    end
                end
            end
        
            % Compute radii (check and correct n_tp if the radii is
            % fixed) - useless at this stage (computed in diffusion1D)
            for iTP = 1:nTP
                Rradi_0{iTP} = ((ME_solid_TP_v0(:,iTP)) ./ ((4/3)*pi.*ME_solid(:,3+nTE+iTP))).^(1/3);
                if diff_nTP(iTP)==0, Rradi_0{iTP}=0*Rradi_0{iTP}; end
                if fix_radi(iTP)~=0
                    Rradi_0{iTP} = fix_radi(iTP)*ones(size(Rradi_0{iTP}));
                    ME_solid(:,3+nTE+iTP) = ME_solid_TP_v0(:,iTP)./((4/3)*pi.*Rradi_0{iTP}.^3);
                end
            end
    
            correct_TE_ind = ones(size(MBC_tp{1}))>0;
            correct_TE_met = -ones(size(MBC_tp{1}))>0;
    
            % This function will not change anything at this stage
            [Gamma_tp, ME_tp, TE_tp, MRadi, MRadi_0, ME_solid] = ...
                diffusion1D(...
                    ME_eq_tp, ME_melt_eq, MBC0_tp, MBC_tp, ...
                    ME_solid, ME_solid0, ME_solid_TP_v0, ME_solid_TP_v0, ME_solid_TP_w0, ... % first timestep
                    ME0_tp, ME_solid_Rho_TP, ME_solid_Rho_TP0, ...
                    Mdiff_tp, correct_TE_ind, correct_TE_met, diff_nTP, diff_nTE, ...
                    fix_radi, nTP, nTE, nc, timestep);
    
            MGamma_e_solid = sum(cat(3, Gamma_tp{:}), 3);
            MGamma_e_fluid = zeros(size(ME_fluid, 1), nTE);
    
            ME_fluid_solid = zeros(size(ME_solid,1),3+nTE);
            ME_fluid_solid(:,1:3) = ME_solid(:,1:3);
    
            for TE_pos = 4:3+nTE
                ME_fluid_solid(:,TE_pos) = interp1(ME_fluid(:,2),ME_fluid(:,TE_pos),ME_solid(:,2),'nearest','extrap');
            end 
    
            ME_solid_TP_v = ME_solid_TP_v0;
            ME_solid_TP_w = ME_solid_TP_w0;
            ME_fluid_TP_v = ME_fluid_TP_v0;
            ME_fluid_TP_w = ME_fluid_TP_w0;
    
            for iTP = 1:nTP
                Error_tp{iTP} = 0 * MBC_tp{iTP};
            end
    
            iter_TE = 0;
    
        else
    
        %% SECOND TIMESTEP ONWARDS
            
        % OPENING
    
            % Recover
    
            ME_solid0 = ME_solid;
            ME_fluid0 = ME_fluid;
            ME_fluid_solid0 = ME_fluid_solid;
    
            ME_solid_Rho_TP0 = ME_solid_Rho_TP;
            ME_Rho_solid0 = ME_Rho_solid;
            ME_Rho_fluid0 = ME_Rho_fluid;
    
            for iTP = 1:nTP
                ME0_tp{iTP} = ME_tp{iTP};
                TE0_tp{iTP} = TE_tp{iTP};
            end
    
            MRadi_0 = MRadi;
    
            ME_solid_TP_v0 = ME_solid_TP_v; % Volume fraction thermo phase
            ME_solid_TP_w0 = ME_solid_TP_w; % Mass fraction thermo phase
            ME_fluid_TP_v0 = ME_fluid_TP_v; % Volume fraction thermo phase
            ME_fluid_TP_w0 = ME_fluid_TP_w; % Mass fraction thermo phase
    
            MEGamma_fluid = interp1(XT(1:nxt+1:end,2),NGamma(1:nxt+1:end,1),ME_fluid0(:,2),'linear','extrap'); MEGamma_fluid = repmat(MEGamma_fluid,1,nTE);
            MEGamma_solid = interp1(XT(1:nxt+1:end,2),NGamma(1:nxt+1:end,1),ME_solid0(:,2),'linear','extrap'); MEGamma_solid = repmat(MEGamma_solid,1,nTE);
    
            MERho_fluid = interp1(XT(1:nxt+1:end,2),NRho0(1:nxt+1:end,2),ME_solid0(:,2),'linear','extrap');
            MERho_solid = interp1(XT(1:nxt+1:end,2),NRho0(1:nxt+1:end,1),ME_solid0(:,2),'linear','extrap');
    
            % Move particles
            NV_mat_TE = NV_mat;
            NV_mat_TE(NV_mat_TE(:,2) > NV_mat_TE(:,4), 4) = NV_mat_TE(NV_mat_TE(:,2) > NV_mat_TE(:,4), 2);
            [ME_solid,ME_fluid] = reactivetransportREE(ME_solid0,ME_fluid0,NV_mat_TE,NDiv,gridxt,gridyt,nxt,0,nyt,0,TT,markmove,timestep,nTP);
            [ME_fluid0_solid,ME_fluid0_fluid] = reactivetransportREE(ME_solid,ME_fluid,[-NV_mat_TE(:,3) -NV_mat_TE(:,4) -NV_mat_TE(:,3) -NV_mat_TE(:,4)],NDiv,gridxt,gridyt,nxt,0,nyt,0,TT,markmove,timestep,nTP);
    
            % Compute TP abundances (mass %)
        
                nSolid = size(ME_solid, 1);
                nFluid = size(ME_fluid, 1);
                ME_solid_TP_w = zeros(nSolid, nTP+1);
                ME_fluid_TP_w = zeros(nFluid, nTP+1);
               
                for iTP = 1:nTP+1
                    ME_solid_TP_w(:,iTP) = interp1(XT(1:nxt+1:end,2),NTp_wt(1:nxt+1:end,iTP),ME_solid(:,2));
                end
                ME_solid_TP_w(ME_solid_TP_w(:,melt_idx)<tol_NTE,melt_idx)=0;
                ME_solid_TP_w = ME_solid_TP_w./repmat(sum(ME_solid_TP_w,2),1,size(ME_solid_TP_w,2));
    
                for iTP = 1:nTP+1
                    ME_fluid_TP_w(:,iTP)  = interp1(XT(1:nxt+1:end,2), NTp_wt(1:nxt+1:end,iTP), ME_fluid(:,2), 'linear', 'extrap');
                end
                ME_fluid_TP_w(ME_fluid_TP_w(:,iTP)<tol_NTE,iTP)=0;
                ME_fluid_TP_w       = ME_fluid_TP_w./repmat(sum(ME_fluid_TP_w,2),1,size(ME_fluid_TP_w,2)); 
    
            % Compute TP abundances (vol %)
        
                nSolid = size(ME_solid, 1);
                nFluid = size(ME_fluid, 1);
                ME_solid_TP_v = zeros(nSolid, nTP+1);
                ME_fluid_TP_v = zeros(nFluid, nTP+1);
            
                for iTP = 1:nTP+1
                    ME_solid_TP_v(:,iTP) = interp1(XT(1:nxt+1:end,2),NTp_vol(1:nxt+1:end,iTP),ME_solid(:,2));
                end
                ME_solid_TP_v(ME_solid_TP_v(:,melt_idx)<tol_NTE,melt_idx)=0;
                ME_solid_TP_v = ME_solid_TP_v./repmat(sum(ME_solid_TP_v,2),1,size(ME_solid_TP_v,2));
    
                ME_solid_TP_v0 = ME_solid_TP_v;
    
                for iTP = 1:nTP+1
                    ME_fluid_TP_v(:,iTP)  = interp1(XT(1:nxt+1:end,2), NTp_vol(1:nxt+1:end,iTP), ME_fluid(:,2), 'linear', 'extrap');
                end
                ME_fluid_TP_v(ME_fluid_TP_v(:,iTP)<tol_NTE,iTP)=0;
                ME_fluid_TP_v       = ME_fluid_TP_v./repmat(sum(ME_fluid_TP_v,2),1,size(ME_fluid_TP_v,2));                 
    
            % Update phi
            FPhi_solid = scatteredInterpolant(XT(:,1), XT(:,2), NPhi(:,1)-NPhi0(:,1), 'nearest', 'nearest'); %Interpolate differences from nodes
            ME_solid(:,3) = ME_solid0(:,3) + FPhi_solid(ME_solid0(:,1), ME_solid0(:,2));
            ME_fluid(:,3) = ME_fluid0(:,3) - FPhi_solid(ME_fluid0(:,1), ME_fluid0(:,2));
    
            ME_fluid(ME_fluid(:,3) < tol_NTE,3) = 0;
    
            % Update densities
            ME_solid_Rho_TP = zeros(nSolid, melt_idx);
            ME_fluid_Rho_TP = zeros(nFluid, melt_idx);
            for iTP = 1:nTP
                ME_solid_Rho_TP(:,iTP) = interp1(XT(1:nxt+1:end,2), NRho_tp(1:nxt+1:end,iTP), ME_solid(:,2));
                ME_fluid_Rho_TP(:,iTP) = interp1(XT(1:nxt+1:end,2), NRho_tp(1:nxt+1:end,iTP), ME_fluid(:,2), 'linear', 'extrap');
            end
            NRho(NRho(:,2)==0,2) = NRho(NRho(:,2)==0,1);
            ME_solid_Rho_TP(:,melt_idx) = interp1(XT(1:nxt+1:end,2), NRho(1:nxt+1:end,2), ME_solid(:,2));
            ME_fluid_Rho_TP(:,melt_idx) = interp1(XT(1:nxt+1:end,2), NRho(1:nxt+1:end,2), ME_fluid(:,2), 'linear', 'extrap');
    
            MERho_fluid = ME_solid_Rho_TP(:,melt_idx);
    
            % Diffusion coefficients & Kd
    
            % P-T in particles;
            ME_T = interp1(XT(1:nxt+1:end,2), NT(1:nxt+1:end,1), ME_solid(:,2)) + constants.TK;
            ME_P = interp1(XT(1:nxt+1:end,2), NP_mat(1:nxt+1:end,1), ME_solid(:,2)) * constants.PGPa;
    
            % Coefficients/diffusivities   
            D_tp = cell(1,nTP);
            E_tp = cell(1,nTP);
            V_tp = cell(1,nTP);
            Kd_tp = cell(1,nTP);
            Mdiff_tp = cell(1,nTP);
    
            for iTP = 1:nTP
                D_tp{iTP}  = repmat(D_coeff(iTP,:), nSolid, 1);
                E_tp{iTP}  = repmat(E_coeff(iTP,:), nSolid, 1);
                V_tp{iTP}  = repmat(V_coeff(iTP,:), nSolid, 1);
                Kd_tp{iTP} = repmat(Kd_coeff(iTP,:), nSolid, 1);
                Mdiff_tp{iTP} = double(D_tp{iTP} .* exp((-E_tp{iTP} + repmat(ME_P,1,nTE) .* V_tp{iTP}) ./ (constants.R .* repmat(ME_T,1,nTE))));
            end
            
            % Calculate Eu diffusivity in Cpx
            if Eu_diff
                iCpx =  find(TP_list=='Cpx');
                if isempty(iCpx); error('Cpx not found in TP_list.'); end
                if isscalar(iCpx)
                    Mdiff_tp{iCpx}  = calcEuCpxdiffusivity(Mdiff_tp{iCpx}, ME_T, ME_P, TE_list, r_Eu, syntheticCpx);
                else
                    for i = 1:length(iCpx)
                        Mdiff_tp{iCpx(i)}  = calcEuCpxdiffusivity(Mdiff_tp{iCpx(i)}, ME_T, ME_P, TE_list, r_Eu, syntheticCpx);
                    end
                end
            end

            % No diffusion if there is not fluid phase nearby
            zeros_aux = zeros(size(ME_solid,1), nTE);
            mask_no_fluid = ME_solid(:,3) > 1 - tol_NTE; % True if (almost) no fluid phase locally
            for iTP = 1:nTP
                Mdiff_tp{iTP}(mask_no_fluid, :) = zeros_aux(mask_no_fluid, :);
            end
    
        % TRACE ELEMENTS
    
            % Boundary conditions
    
            ME_fluid_solid = zeros(size(ME_solid,1), 3+nTE);
            ME_fluid_solid(:,1:2) = ME_solid(:,1:2);         % Copy x, y coordinates
            ME_fluid_solid(:,3)   = 1 - ME_solid(:,3);       % Fluid fraction is 1 - phi_solid
            ME_fluid_solid(ME_fluid_solid(:,3) < tol_NTE, melt_idx) = 0; % If fluid fraction below threshold, set to zero
            
            % Interpolate each TE from fluid points to solid points (by vertical position)
            [~, ia, ~] = unique(ME_fluid(:,2)); % unique y-positions in fluid
            for j = 4 : 3 + nTE
                ME_fluid_solid(:,j) = interp1(ME_fluid(ia,2), ME_fluid(ia, j), ME_solid(:,2), 'nearest', 'extrap');
            end
            
            % (If needed) Advection correction for melt TP in the solid
            ME_RhoTP_solid0_adv = interp1(ME_solid0(:,2), ME_Rho_solid0 .* ME_solid_TP_w0(:, melt_idx), ME_fluid0_solid(:,2), 'linear', 'extrap');
            ME_fluid_solid_adv = ME_fluid_solid;
            
            % Logical index for any nonzero TE in the fluid at each solid point
            ind_fluid = any(ME_fluid_solid(:,4 : 3 + nTE), 2);
    
            % Prepare for next iteration
            for iTP = 1:nTP
                MBC_tp{iTP} = MBC0_tp{iTP}; % Restore BC to previous (or reference) state for all phases
                Kd_tp{iTP} = repmat(Kd_coeff(iTP,:), size(ME_solid,1), 1); % Prepare Kd for all phases
                MBC_tp{iTP}(ind_fluid, :) = ME_fluid_solid(ind_fluid, 4:3+nTE) .* Kd_tp{iTP}(ind_fluid, :); % Update BC where fluid is present, for all phases
            end
    
            % Mass balance constraint: upper and lower bounds for melt's TE concentrations
            for iTP = 1:nTP
                auxE_tp{iTP} = repmat(ME_solid_TP_v(:,iTP) > 0, 1, nTE); % for each mineral phase
            end
            auxE_tp{melt_idx} = repmat(ME_solid_TP_v(:,melt_idx) > 0, 1, nTE); % for melt
    
            % TP volume fractions
            nSolid = size(ME_solid_RhoPhi_tp{1}, 1); 
            nFluid = size(ME_fluid_RhoPhi_tp{1}, 1);
            ME_Rho_solid = zeros(nSolid, 1); 
            ME_Rho_fluid = zeros(nFluid, 1);
            for iTP = 1:nTP+1
                ME_solid_RhoPhi_tp{iTP} = repmat(ME_solid_TP_v(:,iTP) .* ME_solid_Rho_TP(:,iTP), 1, nTE); % Mass TP / vol T
                ME_fluid_RhoPhi_tp{iTP} = repmat(ME_fluid_TP_v(:,iTP) .* ME_fluid_Rho_TP(:,iTP), 1, nTE); % Mass TP / vol T
                ME_Rho_solid = ME_Rho_solid + ME_solid_RhoPhi_tp{iTP}(:,1); % Mass T / vol T
                ME_Rho_fluid = ME_Rho_fluid + ME_fluid_RhoPhi_tp{iTP}(:,1); % Mass T / vol T
            end
        
            % TP mass fractions
            for iTP = 1:nTP+1
                ME_solid_w_tp{iTP} = repmat(ME_solid_TP_w(:,iTP), 1, nTE); % Mass TP / mass T
            end
            ME_solid_w_solid = zeros(nSolid, nTE);
            for iTP = 1:nTP
                ME_solid_w_solid = ME_solid_w_solid + ME_solid_w_tp{iTP}; % % Solid [mass solid / mass T] - normalize to 1 every TP contribution
            end
            ME_solid_w_total = ME_solid_w_solid + ME_solid_w_tp{melt_idx}; % Factor [virtual total mass / total mass]
    
            % Melt TE calculations
    
            denom = zeros(nSolid, nTE);
            for iTP = 1:nTP
                denom = denom + auxE_tp{iTP} .* Kd_tp{iTP} .* ME_solid_w_tp{iTP} ./ ME_solid_w_solid; % [melt] = [solid]/Kd_bulk (normalized to solid)
            end
            Melt_TE = ME_TE ./ denom;
       
            ME_melt_eq = zeros(size(ME_solid,1), nTE);
            meltMask = ME_solid_TP_v(:,melt_idx)>0; % only where melt is present
            ME_melt_eq(meltMask,:) = Melt_TE(meltMask,:); % TE mass / fluid mass
    
            % TE mass / TP mass
            for iTP = 1:nTP
                ME_eq_tp{iTP} = Melt_TE .* auxE_tp{iTP} .* Kd_tp{iTP};
                phaseMask = ME_solid_TP_v(:,iTP)>0;
                ME_tp{iTP} = zeros(nSolid, nTE);
                ME_tp{iTP}(phaseMask,:) = ME_eq_tp{iTP}(phaseMask,:);
            end
    
            correct_TE_met = logical(repmat(ME_solid_TP_v(:,melt_idx)<tol_NTE,1,nTE));
            correct_TE_met_aux = logical(0*correct_TE_met);
    
            % TE profile arrays
            for iTP = 1:nTP
                for TE_pos = 1:nTE
                    TEprofile_tp{iTP,TE_pos} = repmat(ME_eq_tp{iTP}(:,TE_pos)', nc, 1)';
                end
                Mat_aux_tp{iTP} = cell2mat(TEprofile_tp(iTP,:));
                Mat_order = reshape(Mat_aux_tp{iTP}', nc, [])';
                ME_tp{iTP} = mat2cell(Mat_order, nTE*ones(1, nSolid), nc);
            end
    
            % Prepare for iterations
            for iTP = 1:nTP
                Error_tp{iTP} = [];
                MBC_iter_tp{iTP} = [];
            end
    
            MGamma_e_solid_iter = zeros(size(ME_solid,1),nTE);
            MGamma_e_solid_iter_aux = [];
            MGamma_e_solid_max =  zeros(size(ME_solid,1),nTE);
            MGamma_e_solid_min =  zeros(size(ME_solid,1),nTE);
            MGamma_e_solid_max_aux =  zeros(size(ME_solid,1),nTE);
    
            ME_fluid_aux_iter_aux = [];
            ind_gain_loss = logical(zeros(size(ME_solid,1),nTE));
            ind_loss_gain = logical(zeros(size(ME_solid,1),nTE));
    
            TE_Solid_iter = [];
            TE_Fluid_iter = [];
            index_double = [];
            ind_gain = logical(zeros(size(ME_solid,1),nTE));
            ind_loss = logical(zeros(size(ME_solid,1),nTE));
            ind_negative_total = logical(zeros(size(ME_solid,1),nTE));
    
            correct_TE_ind = ones(size(ME_solid,1),nTE)>0;
            correct_TE_ind = correct_TE_ind & correct_TE_met==0;
    
            tol_noChange = tol_TE*tol_NTE/timestep;
    
            for iTP = 1:nTP
                ME_noChange_tp{iTP} = cell2mat(ME_tp{iTP});
                TE_noChange_tp{iTP} = TE_tp{iTP};
            end
    
            ME_solid_noChange = ME_solid(:,4:3+nTE);
            aux_negative_ind = zeros(1, nTP);
    
            % Main TE iterations
    
                iter_TE = 0; exit_TE = 0;
                while exit_TE == 0
                    iter_TE = iter_TE +1;
    
                    correct_TE_met(correct_TE_met_aux) = true;
                    ind_correct_TE_met_aux = correct_TE_met_aux;
    
                    % Compute diffusion
                        % Input:
                            %   -> Boundary conditions t_n and t_{n+1}. 
                            %               E.g. NOli_eq(nodes) and BC_Oli_0(particles)
                            %   -> X,Y coordinates of particles (MR_solid) 
                            %   -> Grain densities (MR_solid)
                            %   -> Proportion (vol%) of TP inside solid (MR_solid_TP_w)
                            %   -> Profiles. E.g. MR_Oli
                            %   -> Diffusion and timestep
                        % Output:
                            %   -> Gradients of every oxide inside TP
                            %   -> New profiles
                            %   -> Boundary conditions for next time step.
    
                            [Gamma_tp, ME_tp, TE_tp, MRadi, MRadi_0, ME_solid] = ...
                                diffusion1D( ...
                                    ME_eq_tp, ME_melt_eq, MBC0_tp, MBC_tp, ...
                                    ME_solid, ME_solid0, ME_solid_TP_v, ME_solid_TP_v0, ME_solid_TP_w, ... % second timestep onwards
                                    ME0_tp, ME_solid_Rho_TP, ME_solid_Rho_TP0, ...
                                    Mdiff_tp, correct_TE_ind, correct_TE_met, diff_nTP, diff_nTE, ...
                                    fix_radi, nTP, nTE, nc, timestep);
    
                    % find where nothing changes
                    if iter_TE == 1
                        Gamma_sum = zeros(size(Gamma_tp{1}));
                        for iTP = 1:nTP
                            Gamma_sum = Gamma_sum + Gamma_tp{iTP};
                        end
                        ind_noChange = abs(Gamma_sum) < tol_noChange;
                        ind_noChange = logical(zeros(size(ind_noChange)));
                    end
                    mat_ind_noChange = reshape(ind_noChange', [], 1);
    
                    % leave unchanged where nothing happens
                    row_cell = nTE * ones(size(cell2mat(ME_tp{1}), 1) / nTE, 1);
                    for iTP = 1:nTP
                        % TE cell arrays
                        mat_ME_tp{iTP} = cell2mat(ME_tp{iTP});
                        mat_ME_tp{iTP}(mat_ind_noChange,:) = ME_noChange_tp{iTP}(mat_ind_noChange,:);    
                        ME_tp{iTP} = mat2cell(mat_ME_tp{iTP}, row_cell, nc);
                        TE_tp{iTP}(ind_noChange) = TE_noChange_tp{iTP}(ind_noChange);
                    end
                
                    ME_solid_noChange_aux = ME_solid(:,4:3+nTE);
                    ME_solid_noChange_aux(ind_noChange) = ME_solid_noChange(ind_noChange);
                    ME_solid(:,4:3+nTE) = ME_solid_noChange_aux;
                                       
                    % Update where something happens
                    if iter_TE == 1
                        % Save current state for all phases (for use as reference in later corrections)
                        mat_ME_tp = cell(1, nTP);
                        for iTP = 1:nTP
                            mat_ME_tp{iTP} = cell2mat(ME_tp{iTP});
                        end
                        Gamma_tp_iter = Gamma_tp;               % Save Gamma for all phases
                        ME_tp_iter    = mat_ME_tp;              % Save ME for all phases (matrix)
                        TE_tp_iter    = TE_tp;                  % Save TE for all phases
                        ME_solid_iter = ME_solid(:,4:3+nTE);    % Save solid TE only
                    else
                        % Correction/update for this iteration
                        for iTP = 1:nTP
                            mat_ME_tp{iTP} = cell2mat(ME_tp{iTP});
                        end
                        ME_solid_aux = ME_solid(:,4:3+nTE);
                    
                        % Correction block for all phases
                        mat_correct_TE_ind = reshape(correct_TE_ind',[],1);
                        Gamma_tp_correct = Gamma_tp_iter;
                        ME_tp_correct    = ME_tp_iter;
                        TE_tp_correct    = TE_tp_iter;
                        ME_solid_correct = ME_solid_iter;
                    
                        for iTP = 1:nTP
                            % Update Gamma
                            Gamma_tp_correct{iTP}(correct_TE_ind) = Gamma_tp{iTP}(correct_TE_ind);
                            % Update ME_tp
                            ME_tp_correct{iTP}(mat_correct_TE_ind,:) = mat_ME_tp{iTP}(mat_correct_TE_ind,:);
                            % Update TE_tp
                            TE_tp_correct{iTP}(correct_TE_ind) = TE_tp{iTP}(correct_TE_ind);
                        end
                        % Update solid TE
                        ME_solid_correct(correct_TE_ind) = ME_solid_aux(correct_TE_ind);
                    
                        % Update for next iteration
                        Gamma_tp_iter = Gamma_tp_correct;
                        ME_tp_iter    = ME_tp_correct;
                        TE_tp_iter    = TE_tp_correct;
                        ME_solid_iter = ME_solid_correct;
                    
                        % Save data back to working variables
                        for iTP = 1:nTP
                            Gamma_tp{iTP} = Gamma_tp_correct{iTP};
                            ME_tp{iTP}    = mat2cell(ME_tp_correct{iTP}, row_cell, nc);
                            TE_tp{iTP}    = TE_tp_correct{iTP};
                        end
                        ME_solid(:,4:3+nTE) = ME_solid_correct;
                    end
    
                    % For flexible gain/loss, loop over TPs and/or operate vectorized.
                    Gamma_sum = zeros(size(Gamma_tp{1}));
                    for iTP = 1:nTP
                        Gamma_sum = Gamma_sum + Gamma_tp{iTP};
                    end
                    MGamma_e_solid = Gamma_sum;  % total solid TE gain/loss for all phases
                    MGamma_e_solid(ind_noChange)=0;
                    MGamma_e_solid_save = MGamma_e_solid;
                
                    % Identify gain/loss
                    ind_gain = MGamma_e_solid - MGamma_e_solid_iter > 0;
                    ind_loss = MGamma_e_solid - MGamma_e_solid_iter < 0;
    
                    if iter_TE > 1
                        % Switch gain/loss flags when sign changes
                        ind_gain_loss(ind_gain_loss == 0 & (ind_gain & ind_loss_iter)) = 1;  % switched from loss to gain
                        ind_loss_gain(ind_loss_gain == 0 & (ind_loss & ind_gain_iter)) = 1;  % switched from gain to loss
                    
                        % Update min/max bounds for solid TE change
                        MGamma_e_solid_min(ind_loss & ind_loss_iter & ind_gain_loss == 0) = MGamma_e_solid(ind_loss & ind_loss_iter & ind_gain_loss == 0);
                        MGamma_e_solid_max(ind_loss & ind_loss_iter) = MGamma_e_solid_iter(ind_loss & ind_loss_iter);
                    
                        MGamma_e_solid_min(ind_gain & ind_loss_iter) = MGamma_e_solid_iter(ind_gain & ind_loss_iter); 
                        MGamma_e_solid_max(ind_loss & ind_gain_iter) = MGamma_e_solid_iter(ind_loss & ind_gain_iter);
                    
                        MGamma_e_solid_min(ind_gain & ind_gain_iter) = MGamma_e_solid_iter(ind_gain & ind_gain_iter);
                        MGamma_e_solid_max(ind_gain & ind_gain_iter & ind_loss_gain == 0) = MGamma_e_solid(ind_gain & ind_gain_iter & ind_loss_gain == 0);
                    
                        % Ensure max bound does not exceed absolute max
                        MGamma_e_solid_max(MGamma_e_solid_max > MGamma_e_solid_max_abs) = MGamma_e_solid_max_abs(MGamma_e_solid_max > MGamma_e_solid_max_abs);
                    
                        % Damped update for melt markers (flexible: use melt index = melt_idx)
                        melt_mask = ME_solid_TP_v(:,melt_idx) > 0;
                        MGamma_e_solid(melt_mask,:) = MGamma_e_solid_min(melt_mask,:) + ...
                            damping0 * (MGamma_e_solid_max(melt_mask,:) - MGamma_e_solid_min(melt_mask,:));
                    end
    
                    % Solid TE gain-loss
                    ind_gain_iter = ind_gain;
                    ind_loss_iter = ind_loss;
    
                    MGamma_e_fluid = MGamma_e_solid;
    
                    % -------------------------------------------------------------------------
                    % PURPOSE: Update trace element concentrations in the fluid phase (ME_fluid_aux)
                    %          using kinetic exchange model with solid phase. Two schemes: 
                    %          exponential solution (for first-order ODEs), or Euler fallback.
                    % -------------------------------------------------------------------------
                    
                    % ---------- a_fluid: Controls exponential decay rate ----------
                    % Legacy/debug code - this sets a_fluid to 0 and does nothing useful
                    a_fluid = 0*MEGamma_solid ./ (repmat(ME_fluid_solid(:,3) .* MERho_fluid,1,nTE)); 
                    a_fluid(isnan(a_fluid)) = 0;
                    
                    % This is the actual (and still zeroed) version of a_fluid
                    a_fluid = 0*MEGamma_solid ./ (repmat(ME_Rho_solid .* ME_solid_TP_w(:,melt_idx),1,nTE));
                    a_fluid(isnan(a_fluid)) = 0;
                    
                    % --> As coded, a_fluid == 0 everywhere, which disables exponential integration
                    
                    % ---------- b_fluid: Source/sink term driving concentration change ----------
                    % First (possibly obsolete) version
                    b_fluid = -damping1 * MGamma_e_fluid ./ (repmat(ME_fluid_solid(:,3) .* MERho_fluid,1,nTE));
                    b_fluid(isnan(b_fluid)) = 0;
                    
                    % Actual version used: damping1 factor scales the exchange term
                    % ME_solid_TP_w(:,melt_idx) is probably melt fraction or mass term
                    b_fluid = -damping1 * MGamma_e_fluid ./ (repmat(ME_Rho_solid .* ME_solid_TP_w(:,melt_idx), 1, nTE));
                    b_fluid(isnan(b_fluid)) = 0;
                    
                    % ---------- c_fluid and c_fluid0: Equilibrium/reference fluid concentrations ----------
                    % Unused: possibly a placeholder for a future implementation
                    c_fluid_aux = ME_fluid_solid0(:,4:3+nTE); % WARNING: unused
                    
                    % Current fluid values at advected locations, rescaled by solid mass/density
                    c_fluid = (repmat(ME_RhoTP_solid0_adv,1,nTE)) .* ME_fluid_solid_adv(:,4:3+nTE) ...
                              ./ (repmat(ME_Rho_solid .* ME_solid_TP_w(:,melt_idx),1,nTE));
                    c_fluid(isnan(c_fluid)) = 0; c_fluid(isinf(c_fluid)) = 0;
                    
                    % Initial fluid values at t=0, rescaled similarly (could be used for time integration consistency)
                    c_fluid0 = (repmat(ME_Rho_solid0 .* ME_solid_TP_w0(:,melt_idx),1,nTE)) .* ME_fluid_solid0(:,4:3+nTE) ...
                               ./ (repmat(ME_Rho_solid .* ME_solid_TP_w(:,melt_idx),1,nTE));
                    c_fluid0(isnan(c_fluid0)) = 0; c_fluid0(isinf(c_fluid0)) = 0;
                    
                    % ---------- Time step ----------
                    dt_fluid = timestep;
                    
                    % ---------- Update ME_fluid_aux: Option 1 — exponential solution ----------
                    % This implements the analytical solution to dC/dt = a*C + b
                    ME_fluid_aux = (b_fluid - exp(-a_fluid*dt_fluid) .* (b_fluid - a_fluid.*c_fluid)) ./ a_fluid;
                    
                    % Euler fallback when a_fluid == 0 to avoid divide-by-zero
                    ME_fluid_aux(a_fluid==0) = c_fluid(a_fluid==0) + dt_fluid * b_fluid(a_fluid==0);
                    
                    % ---------- Final assignment: Euler overwrite (used always) ----------
                    % This line overrides everything above — uses explicit Euler scheme for all values
                    ME_fluid_aux = c_fluid + dt_fluid * b_fluid;
                    
                    % ---------- Compute solid-phase exchange flux for diagnostics or limits ----------
                    % Estimate max absolute gamma (exchange term) for the solid, based on updated c_fluid
                    % Scaled by density and melt fraction, then normalized by dt
                    MGamma_e_solid_max_abs = c_fluid .* (repmat(ME_Rho_solid .* ME_solid_TP_w(:,melt_idx),1,nTE)) / dt_fluid;
    
                    % Clean NaN/Inf
                    ME_fluid_aux(isnan(ME_fluid_aux)) = c_fluid(isnan(ME_fluid_aux));
                    ME_fluid_aux(isinf(ME_fluid_aux)) = c_fluid(isinf(ME_fluid_aux));
                    
                    if iter_TE == 1
                        correct_TE_met_aux = ME_fluid_aux < 0 & c_fluid == 0;
                        correct_TE_met_aux = logical(repmat(any(correct_TE_met_aux,2),1,nTE));
                    end
                    
                    % Clamp negatives where c_fluid is zero
                    ME_fluid_aux(ME_fluid_aux<0 & c_fluid==0) = c_fluid(ME_fluid_aux<0 & c_fluid==0);
                    
                    % Negative correction loop
                    ind_negative = ME_fluid_aux < 0 & correct_TE_ind;
                    ind_negative_total = ind_negative | ind_negative_total;
                    aux_negative_ind(iter_TE) = sum(sum(ind_negative));
    
                    iter_negative = 0;
                    while sum(sum(ind_negative))
                        iter_negative = iter_negative + 1;
                    
                        % Update max for those that are too negative
                        MGamma_e_solid_max_aux(ind_negative) = MGamma_e_solid(ind_negative);
                        MGamma_e_solid_max_aux(ind_negative & (MGamma_e_solid_max_abs < MGamma_e_fluid)) = ...
                            MGamma_e_solid_max_abs(ind_negative & (MGamma_e_solid_max_abs < MGamma_e_fluid));
                        MGamma_e_solid(ind_negative) = MGamma_e_solid_min(ind_negative) + ...
                            damping1/2 * (MGamma_e_solid_max_aux(ind_negative) - MGamma_e_solid_min(ind_negative));
                        MGamma_e_fluid = MGamma_e_solid;
                    
                        % Recompute b_fluid (you may want to keep denominator consistent with above)
                        b_fluid = -damping1 * MGamma_e_fluid ./ (repmat(ME_fluid_solid0(:,3) .* MERho_fluid,1,nTE)); 
                        b_fluid(isnan(b_fluid))=0;  b_fluid(isinf(b_fluid))=0;
                        b_fluid = -damping1 * MGamma_e_fluid ./ (repmat(ME_Rho_solid .* ME_solid_TP_w(:,melt_idx), 1, nTE));
                        b_fluid(isnan(b_fluid)) = 0; b_fluid(isinf(b_fluid)) = 0;
                    
                        ME_fluid_aux(a_fluid==0) = c_fluid(a_fluid==0) + dt_fluid * b_fluid(a_fluid==0);
                    
                        if iter_negative == iter_negative_max
                            ME_fluid_aux(ind_negative) = ME_fluid_aux_iter(ind_negative);
                        end
                    
                        ind_negative = ME_fluid_aux < 0 & correct_TE_ind;
                    end
                    
                    if iter_negative > 0
                        fprintf('Negative fluid composition reached at iteration %d. And exit in %d inner iterations\n', iter_TE, iter_negative);
                    end
    
    
                    if iter_TE == 1
                        % First TE iteration: set min/max for gains/losses
                        MGamma_e_solid_min(ind_loss) = MGamma_e_solid(ind_loss);
                        MGamma_e_solid_max(ind_gain) = MGamma_e_solid(ind_gain);
                    end
                    
                    % Clean up ME_fluid_aux
                    ME_fluid_aux(isnan(ME_fluid_aux)) = c_fluid(isnan(ME_fluid_aux));
                    ME_fluid_aux(isinf(ME_fluid_aux)) = c_fluid(isinf(ME_fluid_aux));
                    ME_fluid_aux(correct_TE_met) = 0;
                    
                    if iter_TE == 1
                        ME_fluid_aux_iter = ME_fluid_aux;
                    else
                        ME_fluid_correct = ME_fluid_aux_iter;
                        ME_fluid_correct(correct_TE_ind) = ME_fluid_aux(correct_TE_ind);
                        ME_fluid_correct(correct_TE_met) = 0;
                        ME_fluid_aux_iter = ME_fluid_correct;
                    end
                    
                    MGamma_e_solid_iter_aux = [MGamma_e_solid_iter_aux, MGamma_e_solid];
                    MGamma_e_solid_iter = MGamma_e_solid;
                    ME_fluid_aux_iter_aux = [ME_fluid_aux_iter_aux, ME_fluid_aux_iter];
                    
                    % Update fluid composition in ME_fluid_solid
                    ME_fluid_solid(:, 4:3+nTE) = ME_fluid_aux_iter;
                    
                    % Prepare BC for next iteration (FLEXIBLE VERSION)
                    mask_fluid = (ME_solid(:,3) < 1-tol_NTE);
                    for iTP = 1:nTP
                        MBC_tp{iTP}(mask_fluid,:) = ME_fluid_solid(mask_fluid,4:3+nTE) .* Kd_tp{iTP}(mask_fluid,:);
                    end
                    
                    % Extract TE compositions per mineral for each TE and phase
                    % Assumes ME_tp is cell array {1,nTP}, each cell is a cell array of size (nGrains,1), each cell nTE*nc
                    for iTP = 1:nTP
                        ME_mat = cell2mat(ME_tp{iTP});
                        for TE_pos = 1:nTE
                            TE = strtrim(TE_list(TE_pos,:)); % e.g. 'La', 'Ce', etc.
                            ME_mat_TE{iTP}.(TE) = ME_mat(TE_pos:nTE:end, :);
                            if TE_pos == 1
                                MBC_aux{iTP} = ME_mat_TE{iTP}.(TE)(:,end);
                            else
                                MBC_aux{iTP} = cat(2, MBC_aux{iTP}, ME_mat_TE{iTP}.(TE)(:,end));
                            end
                        end
                    end
                                               
                    % --- Correct for solid without fluid without TE ---
                    mask_no_fluid = ME_solid(:,3) >= 1-tol_NTE;
                    for iTP = 1:nTP
                        MBC_tp{iTP}(mask_no_fluid,:) = MBC_aux{iTP}(mask_no_fluid,:);
                    end
                    
                    % --- Correct for solid with fluid but without TE ---
                    for iTP = 1:nTP
                        MBC_tp{iTP}(ind_correct_TE_met_aux) = MBC_aux{iTP}(ind_correct_TE_met_aux);
                    end
                    
                    % --- Plotting/prep for diagnostics ---
                    MRadi_plot   = MRadi;   MRadi_plot(MRadi_plot==0) = NaN;
                    MRadi_0_plot = MRadi_0; MRadi_0_plot(MRadi_0_plot==0) = NaN;
                    
                    % --- Flexible sum of all solid TE ---
                    TE_Solid = zeros(size(TE_tp{1}));
                    for iTP = 1:nTP
                        TE_Solid = TE_Solid + TE_tp{iTP};
                    end
                    
                    index_double = [index_double; (iter_TE-1)*size(TE_Solid,1) + find(ind_fluid)];
                    
                    % --- Gather per-phase ME matrices for error computation ---
                    for iTP = 1:nTP
                        mat_ME_tp{iTP} = cell2mat(ME_tp{iTP});
                    end
                    
                    % --- Error computation, phase-by-phase ---
                    Error_iter_tp = cell(1,nTP);                      
                    for iTP = 1:nTP
                        mat_last = reshape(mat_ME_tp{iTP}(:, end), nTE, [])'; % The last column in mat_ME_tp{iTP} is the current value
                        Error_iter{iTP} = (MBC_tp{iTP} - mat_last) ./ MBC_tp{iTP}; % Calculate error (same shape as MBC_tp{iTP}: nSolid x nTE)
                        Error_iter{iTP}(MRadi(:,iTP) == 0, :) = 0; % Zero error where phase is not present  
                        Error_iter{iTP}(correct_TE_met) = 0; % Zero error for masked TEs (e.g. negative or "correct_TE_met")
                        % (Optional: zero error if abs(MBC-mat_last)<tol_TE_abs)
                        % Error_iter{iTP}(abs(MBC_tp{iTP} - mat_last) < tol_TE_abs) = 0;
                    end
    
                    correct_TE_ind_tp = cell(1,nTP);
                    for iTP = 1:nTP
                        nSolid = size(MBC_tp{iTP}, 1);
                        mat_last = reshape(mat_ME_tp{iTP}(:,end), nTE, [])';
                        correct_TE_ind_tp{iTP} = (abs(Error_iter{iTP}) > tol_TE) & (abs(MBC_tp{iTP} - mat_last) > tol_TE_abs);
                    end
                    
                    % Merge all phases to one index
                    correct_TE_ind = false(size(correct_TE_ind_tp{1}));
                    for iTP = 1:nTP
                        correct_TE_ind = correct_TE_ind | correct_TE_ind_tp{iTP};
                    end
                    
                    % --- Apparent Kd (optional, per phase) ---
                    Kd_appa_tp = cell(1,nTP);        
                    for iTP = 1:nTP
                        Kd_appa_tp{iTP} = TE_tp{iTP} ./ ME_fluid_solid(:,4:end);
                        Kd_appa_tp{iTP}(isinf(Kd_appa_tp{iTP})) = NaN;
                        Kd_appa_tp{iTP}(ME_fluid_solid(:,4:end) < tol_TE) = NaN;
                    end
                    
                    % --- Error history for all phases ---
                    for iTP = 1:nTP
                        Error_tp{iTP} = [Error_tp{iTP}  Error_iter_tp{iTP}];
                    end
    
                    % Check convergence
                    if any(any(correct_TE_ind~=0))==0
                        exit_TE = 1;
                        fprintf('=== TIMESTEP %d ===\n', iter)
                        fprintf('Time: %d Ky. Convergence reached at iteration %d.\n',time/1e3/365/24/60/60, iter_TE);
                    end
                    if  iter_TE > max_iter_TE
                        exit_TE = 1;
                        fprintf('=== TIMESTEP %d ===\n', iter)
                        fprintf('Time: %d Ky. Convergence NOT reached. Exit at iteration %d.\n',time/1e3/365/24/60/60, iter_TE);
                    end
                end
    
            for TE_pos = 1:nTE
                [~,ia,~] = unique(ME_fluid_solid(:,2));
                ME_fluid(:,3+TE_pos) = interp1(ME_fluid_solid(ia,2),ME_fluid_aux_iter(ia,TE_pos),ME_fluid(:,2),'nearest','extrap');
            end 
    
        % CLOSING
    
            % Insert solid & fluid
                [ME_solid, ME_fluid_solid, ME_solid_TP_w, ME_solid_TP_v, ME_solid_Rho_TP, ME_tp, TE_tp, MRadi] = ...
                    insertReactiveTE_1D_solid(ME_solid, ME_fluid_solid, ME_solid_TP_w, ME_solid_TP_v, ME_solid_Rho_TP, MRadi, ME_tp, TE_tp, input_TP_v, input_TP_w, gridxt, gridyt);
                [ME_fluid, ME_fluid_TP_w, ME_fluid_TP_v, ME_fluid_Rho_TP] = ...
                    insertReactiveTE_1D_fluid(ME_fluid, ME_fluid_TP_w, ME_fluid_TP_v, ME_fluid_Rho_TP, gridxt, gridyt, TE_input_liquid, input_TP_v);
    
            for iTP = 1:nTP
                mat_ME_tp{iTP} = cell2mat(ME_tp{iTP}); % (nTE*ngrain, nc)
                % Take the final col (end) for last time step, then reshape to (nSolid, nTE)
                MBC0_tp{iTP} = reshape(mat_ME_tp{iTP}(:,end), nTE, [])';
            end
    
        end
    
        %% COMPUTE BULK WITH NEW PARTICLES
        
        % Number of solid markers
        nSolid = size(ME_solid,1);
        nFluid = size(ME_fluid,1);
        
        % mass TP / mass T
        % Each phase (including melt) gets a weight fraction matrix (nSolid x nTE)
        for iTP = 1:nTP+1
            ME_solid_w_tp{iTP} = repmat(ME_solid_TP_w(:,iTP),1,nTE);
        end
        
        % solid [mass solid / mass T]
        % Sum weight contributions of all solid phases (not including melt)
        ME_solid_w_solid = zeros(nSolid,nTE);
        for iTP = 1:nTP
            ME_solid_w_solid = ME_solid_w_solid + ME_solid_w_tp{iTP};
        end
        
        % factor [virtual total mass / total mass]
        ME_solid_w_total = ME_solid_w_solid + ME_solid_w_tp{melt_idx}; % includes melt
        
        % --- bulk TE [mass TE/ mass T]
        ME_TE = zeros(nSolid, nTE);
        for iTP = 1:nTP
            ME_TE = ME_TE + ME_solid_w_tp{iTP} .* TE_tp{iTP};
        end
        ME_TE = (ME_TE + ME_solid_w_tp{melt_idx} .* ME_fluid_solid(:,4:3+nTE)) ./ ME_solid_w_total;
        
        % mass TP / vol T
        for iTP = 1:nTP+1
            ME_solid_RhoPhi_tp{iTP} = repmat(ME_solid_TP_v(:,iTP) .* ME_solid_Rho_TP(:,iTP), 1, nTE);
        end
        
        % mass T / vol T
        ME_Rho_solid = zeros(nSolid,1);
        for iTP = 1:nTP+1
            ME_Rho_solid = ME_Rho_solid + ME_solid_RhoPhi_tp{iTP}(:,1);
        end
        
        % mass TP / vol T
        for iTP = 1:nTP+1
            ME_fluid_RhoPhi_tp{iTP} = repmat(ME_fluid_TP_v(:,iTP) .* ME_fluid_Rho_TP(:,iTP), 1, nTE);
        end
        
        % mass T / vol T
        ME_Rho_fluid = zeros(nFluid,1);
        for iTP = 1:nTP+1
            ME_Rho_fluid = ME_Rho_fluid + ME_fluid_RhoPhi_tp{iTP}(:,1);
        end
        
        %% SAVE & UPDATE TIME
        
        time = time + timestep;
        timeMy = time/1e6/365/24/60/60;
        timeMy_mat = [timeMy_mat timeMy];
        
        % Sorting
        ME_solid_sorted = ME_solid(:,2:end);
        ME_fluid_sorted = ME_fluid_solid(:,2:end);
        TE_solid = ME_solid_sorted(:,3:2+nTE);
        TE_fluid = ME_fluid_sorted(:,3:2+nTE);
        
        % Saving
        iter_TE_save = single(iter_TE);
        NT_save = single(NT);
        NPhi_save = single(NPhi);
        NRho_save = single(NRho);
        NRho_tp_save = single(NRho_tp);
        ME_solid_TP_v_save = single(ME_solid_TP_v);
        ME_solid_TP_w_save = single(ME_solid_TP_w);
        ME_solid_Rho_TP_save = single(ME_solid_Rho_TP);
        ME_TE_save = single(ME_TE);
        NV_mat_save = single(NV_mat);
        XT_save = single(XT);
        ME_solid_save = single(ME_solid);
        ME_fluid_save = single(ME_fluid);
        ME_fluid_solid_save = single(ME_fluid_solid);
        MRadi_save = single(MRadi);
        timestep_save = single(timestep);
        time_save = single(time);
        timeMy_save = single(timeMy);
        timeMy_mat_save = single(timeMy_mat);
        
        iter_save = single(iter);
        gridxt_save = single(gridxt);
        gridyt_save = single(gridyt);
        
        % Save per-phase arrays if needed (convert cell to matrix)
        ME_save_tp = cell(1, nTP);
        for iTP = 1:nTP
            ME_save_tp{iTP} = single(cell2mat(ME_tp{iTP}));
        end
        
        aux_mesh = zeros(1, nc);
        for index_mesh = 1:nc
            aux_mesh(index_mesh) = (1-(1/index_mesh)^exp_factor)/(1-(1/nc)^exp_factor);
        end
        for iTP = 1:nTP
            xmesh_tp{iTP} = MRadi(1,iTP)*aux_mesh;
        end
    
        if ismember(iter,save_list)
            baseFileName = sprintf('savetime_%d', save_iter);
            fullFileName = fullfile('output', baseFileName);
            save(fullFileName,'save_list','xmesh_tp','ME_solid_TP_v_save','NV_mat_save','XT_save','ME_solid_save','ME_fluid_save','ME_fluid_solid_save','ME_save_tp', 'MRadi_save','timestep_save','time_save','timeMy_save','iter_save','gridxt_save','gridyt_save','timeMy_mat_save');
            save_iter = save_iter+1;
        end
        
        ME_solid_sorted = ME_solid_save(:,2:end);
        ME_fluid_sorted = ME_fluid_solid_save(:,2:end);
        
        disp('Saving completed')
        iter = iter + 1;
    
    end