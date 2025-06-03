function [ME,ME_fluid_TP_w,ME_fluid_TP_v,ME_fluid_Rho_TP,input_aux] = insertReactiveTE_1D_fluid(ME,ME_fluid_TP_w,ME_fluid_TP_v,ME_fluid_Rho_TP,gridx,gridy,input_REE,input_TP_v)

input_aux = 0;
% Remove particle
aux_correct = ME(:,1)>gridx(end) | ME(:,1)<gridx(1); 
ME(aux_correct,:) = [];
ME_fluid_TP_w(aux_correct,:) = [];
ME_fluid_TP_v(aux_correct,:) = [];
ME_fluid_Rho_TP(aux_correct,:) = [];

aux_correct = ME(:,2)>gridy(end) | ME(:,2)<gridy(1); 
ME(aux_correct,:) = [];
ME_fluid_TP_w(aux_correct,:) = [];
ME_fluid_TP_v(aux_correct,:) = [];
ME_fluid_Rho_TP(aux_correct,:) = [];
% 
% xsize = gridx(end)-gridx(1);
% ysize = gridy(end)-gridy(1);
% ystp = gridy(2)-gridy(1);
% xstp = gridx(2)-gridx(1);
% 
% mxnum   =   1;         % total number of markers in horizontal direction
% mynum   =   1;         % total number of markers in vertical direction
% mxstep  =   xsize/mxnum;        % step between markers in horizontal direction
% mystep  =   ystp/5;        % step between markers in vertical direction
% 
% if min(ME(:,2))>ystp/5
%     stringx = mxstep/2:mxstep:mxnum*mxstep-mxstep/2;
%     MX = repmat(stringx,mynum,1);
%     MX_rand = (rand(size(MX))-0.5)*mxstep;
%     stringy = mystep/2:mystep:mynum*mystep-mystep/2;
%     stringy = 0;
%     MY = repmat(stringy',1,mxnum);
%     MY_rand = (rand(size(MY))-0.5)*mystep;
%     MX_aux = [MX(:)+0*MX_rand(:) MY(:)+0*MY_rand(:)];
%      % add new particles
%     ME = [MX_aux input_TP_v(end)*ones(length(stringx),1) repmat(input_REE,length(stringx),1) repmat(ME(1,4+length(input_REE):end),length(stringx),1); ME ];
%     input_aux = 1;
% end

%%
% Updating temperature for markers
MS_aux = ME;
X_MS = MS_aux(:,1);
Y_MS = MS_aux(:,2);
nxt = length(gridx)-1;
nyt = length(gridy)-1;
xsize = gridx(end)-gridx(1);
ysize = gridy(end)-gridy(1);
ystp = gridy(2)-gridy(1);
xstp = gridx(2)-gridx(1);

% % xn = zeros(size(X_MS));
% % xn(X_MS<gridx(nxt+1))=double(int16(X_MS(X_MS<gridx(nxt+1))./xstp-0.5))+1;
% yn=double(int16(Y_MS./ystp-0.5))+1;
% yn(Y_MS<gridy(nyt+1))=double(int16(Y_MS(Y_MS<gridy(nyt+1))./ystp-0.5))+1;
yn = max(cumsum(Y_MS./gridy>1,2),[],2);

% if (xn<1)
%     xn  =   1;
% end
% if (xn>(nxt))
%     xn  =   (nxt);
% end
if (yn<1)
    yn  =   1;
end
if (yn>(nyt))
    yn  =   (nyt);
end


% remove_index = [];
% MS_add = [];
N_particles = accumarray(yn,ones(size(yn)));
correct_elem = (N_particles<10);
index_ielem = find(correct_elem==1);

if sum(correct_elem)>0; input_aux = 1; end;

mxnum   =   1;         % total number of markers in horizontal direction
mynum_elem = 5;
mynum_elem = 10;
mynum   =   (length(gridy)-1)*mynum_elem;         % total number of markers in vertical direction
mynum   =   (length(gridy)-1);         % total number of markers in vertical direction
mystep  =   ysize/(mynum_elem*mynum);        % step between markers in vertical direction
% mystep  =   gridy(2:end)-gridy(1:end-1);        % step between markers in vertical direction
stringx = 0;
MX = repmat(stringx,mynum,1);
MX_rand = (rand(size(MX))-0.5)*0;
% stringy = mystep/2:mystep:mynum*mystep-mystep/2; 
stringy = gridy(2:end)-(gridy(2:end)-gridy(1:end-1));
MY = repmat(stringy',1,mxnum);
MY_rand = (rand(size(MY))-0.5).*mystep'/5;
MY_rand = (rand(size(MY))-0.5).*mystep;
MY_rand = -mystep/2.*ones(size(MY));
MS_aux = [MX(:)+0*MX_rand(:) MY(:)+MY_rand(:)];
yn_aux = reshape(repmat(1:length(gridy)-1,mynum_elem,1),(length(gridy)-1)*mynum_elem,[]);
yn_aux = reshape(repmat(1:length(gridy)-1,1,1),(length(gridy)-1)*1,[]);

% % find closest particle to copy in the desired value of stringy
% insert_coord = MS_aux(correct_elem);
% insert_coord_mat = repmat(insert_coord,length(Y_MS),1);
% Y_MS_mat = repmat(Y_MS,1,length(insert_coord));
% diff = insert_coord_mat-Y_MS_mat;
% [~, ind_add] = min(abs(diff));


[ind_add, ~] = ismember(yn_aux,index_ielem);
[ind_remove, ~] = ismember(yn,index_ielem); ind_remove = [];
nREE = length(input_REE);
MS_add = [MS_aux(ind_add,:) zeros(size(MS_aux(ind_add,:),1),size(ME,2)-2)];
ME_fluid_TP_w_add = zeros(size(MS_aux(ind_add,:),1),size(ME_fluid_TP_w,2));
ME_fluid_TP_v_add = zeros(size(MS_aux(ind_add,:),1),size(ME_fluid_TP_v,2));
ME_fluid_Rho_TP_add = zeros(size(MS_aux(ind_add,:),1),size(ME_fluid_Rho_TP,2));
for index = 1:size(ME,2)-2
    [~,ia,~] = unique(Y_MS);
    MS_add(:,2+index) = interp1(Y_MS(ia),ME(ia,2+index),MS_add(:,2),'nearest','extrap');
end
MS_add(:,4:3+nREE) = repmat(input_REE,size(MS_add,1),1);

for index = 1:size(ME_fluid_TP_w,2)
    [~,ia,~] = unique(Y_MS);
    ME_fluid_TP_w_add(:,index) = interp1(Y_MS(ia),ME_fluid_TP_w(ia,index),MS_add(:,2),'nearest','extrap');
    ME_fluid_TP_v_add(:,index) = interp1(Y_MS(ia),ME_fluid_TP_v(ia,index),MS_add(:,2),'nearest','extrap');
    ME_fluid_Rho_TP_add(:,index) = interp1(Y_MS(ia),ME_fluid_Rho_TP(ia,index),MS_add(:,2),'nearest','extrap');
end


ME(ind_remove,:) = [];
ME_fluid_TP_w(ind_remove,:) = [];
ME_fluid_TP_v(ind_remove,:) = [];
ME_fluid_Rho_TP(ind_remove,:) = [];
ME = [ME; MS_add];
ME_fluid_TP_w = [ME_fluid_TP_w; ME_fluid_TP_w_add];
ME_fluid_TP_v = [ME_fluid_TP_v; ME_fluid_TP_v_add];
ME_fluid_Rho_TP = [ME_fluid_Rho_TP; ME_fluid_Rho_TP_add];

% Remove the repeated values
ME_preDouble    =  ME(:,2);
[~,ia,~] = unique(ME_preDouble);
ME = ME(ia,:);
ME_fluid_TP_w = ME_fluid_TP_w(ia,:);
ME_fluid_TP_v = ME_fluid_TP_v(ia,:);
ME_fluid_Rho_TP = ME_fluid_Rho_TP(ia,:);

% get the indexes in order
ME_preSorted    =  ME(:,2);
[~, ind_orden] = sort(ME_preSorted,1); 

ME = ME(ind_orden,:);
ME_fluid_TP_w = ME_fluid_TP_w(ind_orden,:);
ME_fluid_TP_v = ME_fluid_TP_v(ind_orden,:);
ME_fluid_Rho_TP = ME_fluid_Rho_TP(ind_orden,:);

% ME = sort([ME(:,2) ME(:,1) ME(:,3:end)],1); 
% ME = [ME(:,2) ME(:,1) ME(:,3:end)];

% 
%     
% for index = 1:sum(correct_elem)
% %     Xe = XT(T(index_ielem(index),:),:);
% %     xsize = max(Xe(:,1))-min(Xe(:,1));
% %     ysize = max(Xe(:,2))-min(Xe(:,2));
%     ysize = gridy(yn(index)+1)-gridy(yn(index));
% %     mxstep  =   xsize/mxnum;        % step between markers in horizontal direction
%     mxstep  = 0;
%     mystep  =   ysize/mynum;        % step between markers in vertical direction
% %     stringx = min(Xe(:,1))+mxstep/2:mxstep:min(Xe(:,1))+mxnum*mxstep-mxstep/2;
%     stringx = 0;
%     MX = repmat(stringx,mynum,1);
%     MX_rand = (rand(size(MX))-0.5)*mxstep;
%     stringy = gridy(yn(index))+mystep/2:mystep:gridy(yn(index))+mynum*mystep-mystep/2;
%     MY = repmat(stringy',1,mxnum);
%     MY_rand = (rand(size(MY))-0.5)*mystep;
%     MX_aux = [MX(:)+MX_rand(:) MY(:)+0*MY_rand(:)];
%     
%     % SiO2[wt%] - Al2O3[wt%] - FeO[wt%] - MgO[wt%] - CaO[wt%]
%     interp1(ME_fluid_solid(:,2),ME_fluid_aux_iter(:,index),ME(:,2),'nearest','extrap');
%     FR_SiO2 = scatteredInterpolant(Xe(:,1),Xe(:,2),NR(T(index_ielem(index),:),1));
%     MR_SiO2_aux = FR_SiO2(MX_aux(:,1),MX_aux(:,2));
%     
%     FR_Al2O3 = scatteredInterpolant(Xe(:,1),Xe(:,2),NR(T(index_ielem(index),:),2));
%     MR_Al2O3_aux = FR_Al2O3(MX_aux(:,1),MX_aux(:,2));
%     
%     FR_FeO = scatteredInterpolant(Xe(:,1),Xe(:,2),NR(T(index_ielem(index),:),3));
%     MR_FeO_aux = FR_FeO(MX_aux(:,1),MX_aux(:,2));
%     
%     FR_MgO = scatteredInterpolant(Xe(:,1),Xe(:,2),NR(T(index_ielem(index),:),4));
%     MR_MgO_aux = FR_MgO(MX_aux(:,1),MX_aux(:,2));
%     
%     FR_CaO = scatteredInterpolant(Xe(:,1),Xe(:,2),NR(T(index_ielem(index),:),5));
%     MR_CaO_aux = FR_CaO(MX_aux(:,1),MX_aux(:,2));
%     
%     FR_Na2O = scatteredInterpolant(Xe(:,1),Xe(:,2),NR(T(index_ielem(index),:),6));
%     MR_Na2O_aux = FR_Na2O(MX_aux(:,1),MX_aux(:,2));
%     
%     FR_Cr2O3 = scatteredInterpolant(Xe(:,1),Xe(:,2),NR(T(index_ielem(index),:),7));
%     MR_Cr2O3_aux = FR_Cr2O3(MX_aux(:,1),MX_aux(:,2));
%     
%     F_Phi = scatteredInterpolant(Xe(:,1),Xe(:,2),NPhi(T(index_ielem(index),:)));
%     MR_Phi_aux = F_Phi(MX_aux(:,1),MX_aux(:,2));
%     
%     % remove particles
%     remove_index = [remove_index; find(ielem==index_ielem(index))];
%     %     MS(ielem==index_ielem(index),:) = [];
%     
%     % add new particles
%     MS_add = [MS_add; MX_aux MR_Phi_aux MR_SiO2_aux MR_Al2O3_aux MR_FeO_aux MR_MgO_aux MR_CaO_aux MR_Na2O_aux MR_Cr2O3_aux zeros(length(MR_CaO_aux),nTP)];
% end
% 
% MR(remove_index,:) = [];
% MR = [MR; MS_add];


