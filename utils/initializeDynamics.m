function [NR_solid,NR_fluid,NGamma,NRho0,NRho,NRho_tp,NRho_tp0,NV_mat0,NV_mat,NP_mat,NDiv0,NDiv,NPhi0,NPhi,NTp_wt0,NTp_wt,NTp_vol0,NTp_vol,NT,timestep,time,label] = ...
    initializeDynamics(fullFileName,nxt,nyt)

load(fullFileName);

NR_solid = NR_solid_save(1:(nxt+1)*(nyt+1),:);
NR_fluid = NR_fluid_save(1:(nxt+1)*(nyt+1),:);
NGamma = NGamma_save(1:(nxt+1)*(nyt+1),:);
NRho0 = NRho0_save(1:(nxt+1)*(nyt+1),:);
NRho = NRho_save(1:(nxt+1)*(nyt+1),:);
NRho_tp = NRho_tp_save(1:(nxt+1)*(nyt+1),:);
NRho_tp0 = NRho_tp_save(1:(nxt+1)*(nyt+1),:);
NV_mat0 = NV_mat0_save(1:(nxt+1)*(nyt+1),:);
NV_mat = NV_mat_save(1:(nxt+1)*(nyt+1),:);
NP_mat = NP_mat_save(1:(nxt+1)*(nyt+1),:);
NDiv0 = NDiv0_save(1:(nxt+1)*(nyt+1),:);
NDiv = NDiv_save(1:(nxt+1)*(nyt+1),:);
NPhi0 = NPhi0_save(1:(nxt+1)*(nyt+1),:);
NPhi = NPhi_save(1:(nxt+1)*(nyt+1),:);
NTp_wt0 = NTp_wt0_save(1:(nxt+1)*(nyt+1),:);
NTp_wt = NTp_wt_save(1:(nxt+1)*(nyt+1),:);
NTp_vol0 = NTp_vol0_save(1:(nxt+1)*(nyt+1),:);
NTp_vol = NTp_vol_save(1:(nxt+1)*(nyt+1),:);
NT = NT_save(1:(nxt+1)*(nyt+1),:);
%MT = MT_save(1:(nxt+1)*(nyt+1),:);

% Delete if there are non-matching zeros
NPhi0(NPhi0(:,2)==0 | NTp_wt0(:,7)==0 | NTp_vol0(:,7)==0,2) = 0;
NTp_wt0(NPhi0(:,2)==0 | NTp_wt0(:,7)==0 | NTp_vol0(:,7)==0,7) = 0;
NTp_vol0(NPhi0(:,2)==0 | NTp_wt0(:,7)==0 | NTp_vol0(:,7)==0,7) = 0;
NPhi(NPhi(:,2)==0 | NTp_wt(:,7)==0 | NTp_vol(:,7)==0,2) = 0;
NTp_wt(NPhi(:,2)==0 | NTp_wt(:,7)==0 | NTp_vol(:,7)==0,7) = 0;
NTp_vol(NPhi(:,2)==0 | NTp_wt(:,7)==0 | NTp_vol(:,7)==0,7) = 0;

for index = 1:size(NTp_wt,2)-1
NTp_wt0(NTp_wt0(:,index)==0 | NTp_vol0(:,index)==0,index) = 0;
NTp_vol0(NTp_wt0(:,index)==0 | NTp_vol0(:,index)==0,index) = 0;
NTp_wt(NTp_wt(:,index)==0 | NTp_vol(:,index)==0,index) = 0;
NTp_vol(NTp_wt(:,index)==0 | NTp_vol(:,index)==0,index) = 0;
end

% timestep = timestep_save/
% time = time_save;
% label = [0 0 0 0 0 1];


