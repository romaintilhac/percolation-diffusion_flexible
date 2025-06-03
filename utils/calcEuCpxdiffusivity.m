function Mdiff_Cpx = calcEuCpxdiffusivity(Mdiff_Cpx, ME_T, ME_P, TE_list, r_Eu, method)
% calculate Eu diffusivities in Cpx based on Sr diffusivities from Sneeringer et al. (1984)

R = constants.R;

Eu_pos = find(strcmp(TE_list, 'Eu'));
if isempty(Eu_pos)
    error('Eu not found in TE_list.');
end
nTE=length(TE_list);

switch method
    case 1 %synthetic diopside
        D_Cpx_Eu2_ref = 1200*1e-4; % [m^2/s]
        E_Cpx_Eu2_ref = 122*4.184*1e3; % [J/mol]
    case 0 %natural diopside
        D_Cpx_Eu2_ref = 54*1e-4; % [m^2/s]
        E_Cpx_Eu2_ref = 97*4.184*1e3; % [J/mol]
end

D_Cpx_Eu2 = repmat(D_Cpx_Eu2_ref,size(ME_T,1),1); 
E_Cpx_Eu2 = repmat(E_Cpx_Eu2_ref,size(ME_T,1),1);
Sn84a=repmat(2.94,size(ME_T,1),1); 
Sn84b=repmat(4640,size(ME_T,1),1);
Mdiff_Cpx_Eu2   = double(D_Cpx_Eu2.*exp(-E_Cpx_Eu2./(R.*ME_T)).*exp(-(ME_P.*1e-9)./((Sn84a.*ME_T-Sn84b).*R))); % [m^2/s]   
Mdiff_Cpx_Eu = Mdiff_Cpx_Eu2.*r_Eu + Mdiff_Cpx(:,Eu_pos).*(1-r_Eu); 
    
Mdiff_Cpx = [Mdiff_Cpx(:,1:Eu_pos-1) Mdiff_Cpx_Eu Mdiff_Cpx(:,Eu_pos+1:nTE)];