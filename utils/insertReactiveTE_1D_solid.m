function [ME_solid, ME_fluid_solid, ME_solid_TP_w, ME_solid_TP_v, ME_solid_Rho_TP, ME_tp, TE_tp, MRadi, input_aux] = ...
    insertReactiveTE_1D_solid(ME_solid, ME_fluid_solid, ME_solid_TP_w, ME_solid_TP_v, ME_solid_Rho_TP, MRadi, ME_tp, TE_tp, input_TP_v, input_TP_w, gridx, gridy)

nTP = size(ME_tp, 2);

input_aux = 0;
% Remove particle
aux_correct = ME_solid(:,1)>gridx(end) | ME_solid(:,1)<gridx(1);
ME_solid(aux_correct,:) = [];
ME_fluid_solid(aux_correct,:) = [];
ME_solid_TP_w(aux_correct,:) = [];
ME_solid_TP_v(aux_correct,:) = [];
ME_solid_Rho_TP(aux_correct,:) = [];
MRadi(aux_correct,:) = [];

for iTP = 1:nTP
    ME_tp{iTP}(aux_correct) = [];
    TE_tp{iTP}(aux_correct,:) = [];
end

aux_correct = ME_solid(:,2)>gridy(end) | ME_solid(:,2)<gridy(1);
ME_solid(aux_correct,:) = [];
ME_fluid_solid(aux_correct,:) = [];
ME_solid_TP_w(aux_correct,:) = [];
ME_solid_TP_v(aux_correct,:) = [];
ME_solid_Rho_TP(aux_correct,:) = [];
MRadi(aux_correct,:) = [];

for iTP = 1:nTP
    ME_tp{iTP}(aux_correct) = [];
    TE_tp{iTP}(aux_correct,:) = [];
end

ystp = gridy(2)-gridy(1);

% Updating temperature for markers
MS_aux = ME_solid;
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
correct_elem = (N_particles<8);

if sum(correct_elem)>0; input_aux = 1; end

mxnum   =   1;         % total number of markers in horizontal direction
mynum_elem = 1;
mynum   =   (length(gridy)-1)*mynum_elem;         % total number of markers in vertical direction
mystep  =   ysize/mynum;        % step between markers in vertical direction
mystep  =   gridy(2:end)-gridy(1:end-1);        % step between markers in vertical direction
stringx = 0;
MX = repmat(stringx,mynum,1);
MX_rand = (rand(size(MX))-0.5)*0;
% stringy = mystep/5:mystep:mynum*mystep-mystep/2;
stringy = gridy(1:end-1)+(gridy(2:end)-gridy(1:end-1))/5;
MY = repmat(stringy',1,mxnum);
MY_rand = (rand(size(MY))-0.5).*mystep'/5;
MS_aux = [MX(:)+0*MX_rand(:) MY(:)+0*MY_rand(:)];

% find closest particle to copy in the desired value of stringy
insert_coord = stringy(correct_elem);
insert_coord_mat = repmat(insert_coord,length(Y_MS),1);
Y_MS_mat = repmat(Y_MS,1,length(insert_coord));
diff = insert_coord_mat-Y_MS_mat;
[~, ind_add] = min(abs(diff));

ME_solid = [ME_solid;...
      ME_solid(ind_add,1) insert_coord' ME_solid(ind_add,3:end)];
ME_fluid_solid = [ME_fluid_solid;...
                      ME_fluid_solid(ind_add,1) insert_coord' ME_fluid_solid(ind_add,3:end)];
ME_solid_TP_w = [ME_solid_TP_w;
                     ME_solid_TP_w(ind_add,:)];
ME_solid_TP_v = [ME_solid_TP_v;...
                     ME_solid_TP_v(ind_add,:)];
ME_solid_Rho_TP = [ME_solid_Rho_TP;...
                       ME_solid_Rho_TP(ind_add,:)];
MRadi = [MRadi;...
         MRadi(ind_add,:)];

for iTP = 1:nTP
    ME_tp{iTP} = [ME_tp{iTP};...
                 ME_tp{iTP}(ind_add,:)];
    TE_tp{iTP} = [TE_tp{iTP};...
                 TE_tp{iTP}(ind_add,:)];
end

% Remove the repeated values
ME_preDouble    =  ME_solid(:,2);
[~,ia,~] = unique(ME_preDouble);
ME_solid = ME_solid(ia,:);
ME_fluid_solid  = ME_fluid_solid(ia,:);
ME_solid_TP_w   = ME_solid_TP_w(ia,:);
ME_solid_TP_v   = ME_solid_TP_v(ia,:);
ME_solid_Rho_TP = ME_solid_Rho_TP(ia,:);
MRadi   = MRadi(ia,:);

for iTP = 1:nTP
    ME_tp{iTP} = ME_tp{iTP}(ia,:);
    TE_tp{iTP} = TE_tp{iTP}(ia,:);
end

% get the indexes in order
ME_preSorted    =  ME_solid(:,2);
[~, ind_orden] = sort(ME_preSorted,1); 
ME_solid = ME_solid(ind_orden,:); 
ME_fluid_solid  = ME_fluid_solid(ind_orden,:);
ME_solid_TP_w   = ME_solid_TP_w(ind_orden,:);
ME_solid_TP_v   = ME_solid_TP_v(ind_orden,:);
ME_solid_Rho_TP = ME_solid_Rho_TP(ind_orden,:);
MRadi   = MRadi(ind_orden,:);

for iTP = 1:nTP
    ME_tp{iTP} = ME_tp{iTP}(ind_orden,:);
    TE_tp{iTP} = TE_tp{iTP}(ind_orden,:);
end
