function [MR_solid,MR_fluid] = reactivetransportREE(MR_solid,MR_fluid,NV,NDiv,gridx,gridy,nxt,nxt_aux,nyt,nyt_aux,T,markmove,dt,nTP)
%
% Moving Markers by velocity field
% Runge-Kutta 1,2,3,4
%
%%

xstp = gridx(2)-gridx(1);
ystp = gridy(2)-gridy(1);

MR_C_solid_aux = MR_solid;
MR_C_fluid_aux = MR_fluid;
NV_aux   = NV;

% add a little bit of porosity to the fluid and solid phases
% MR_solid(MR_solid(:,3)<=1E-4,3)=1E-4;
% MR_fluid(MR_fluid(:,3)<=1E-4,3)=1E-4;

if(markmove>0)
    
    MC0 = MR_C_solid_aux(:,3:5); % expandir
    NVx_solid = NV_aux(:,1);
    NVy_solid = NV_aux(:,2);
    NVx_fluid = NV_aux(:,3);
    NVy_fluid = NV_aux(:,4);
    NDiv_solid = NDiv(:,1);
    X_MS_solid = MR_C_solid_aux(:,1);
    Y_MS_solid = MR_C_solid_aux(:,2);
    X_MS_fluid = MR_C_fluid_aux(:,1);
    Y_MS_fluid = MR_C_fluid_aux(:,2);
    cur_X_solid = X_MS_solid;
    cur_Y_solid = Y_MS_solid;
    cur_X_fluid = X_MS_fluid;
    cur_Y_fluid = Y_MS_fluid;
    
    % Create arrays for velocity of markers
    MVx_solid = zeros(size(X_MS_solid,1),markmove);
    MVy_solid = zeros(size(X_MS_solid,1),markmove);
    MVx_fluid = zeros(size(X_MS_fluid,1),markmove);
    MVy_fluid = zeros(size(X_MS_fluid,1),markmove);
    MDiv_solid = zeros(size(X_MS_solid,1),markmove);
    
    
    for rk=1:1:markmove
%         xn_solid = zeros(size(cur_X_solid));
%         xn_solid(X_MS_solid<gridx(nxt+1))=double(int16(X_MS_solid(X_MS_solid<gridx(nxt+1))./xstp-0.5))+1;
        xn_solid = max(cumsum(X_MS_solid./gridx>1,2),[],2);
        if nxt_aux > 0
            xstp_aux = gridx(end)-gridx(end-1);
            xn_solid(X_MS_solid>gridx(nxt+1))=double(int16((X_MS_solid(X_MS_solid>gridx(nxt+1))-gridx(nxt+1))./xstp_aux-0.5))+1+nxt;
        end
%         yn_solid=double(int16(Y_MS_solid./ystp-0.5))+1;
%         yn_solid(Y_MS_solid<gridy(nyt+1))=double(int16(Y_MS_solid(Y_MS_solid<gridy(nyt+1))./ystp-0.5))+1;
        yn_solid = max(cumsum(Y_MS_solid./gridy>1,2),[],2);
        if nyt_aux > 0
            ystp_aux = gridy(end)-gridy(end-1);
            yn_solid(Y_MS_solid>gridy(nyt+1))=double(int16((Y_MS_solid(Y_MS_solid>gridy(nyt+1))-gridy(nyt+1))./ystp_aux-0.5))+1+nyt;
        end
        if any(xn_solid<1)
            xn_solid(xn_solid<1)  =   1;
        end
        if any(xn_solid>(nxt+nxt_aux))
            xn_solid(xn_solid>(nxt+nxt_aux))  =   (nxt+nxt_aux);
        end
        if any(yn_solid<1)
            yn_solid(yn_solid<1)  =   1;
        end
        if any(yn_solid>(nyt+nyt_aux))
            yn_solid(yn_solid>(nyt+nyt_aux))  =   (nyt+nyt_aux);
        end
        % Element en el que esta ielem = (xnum-1)*(yn-1)+xn
        ielem_solid = (nxt+nxt_aux)*(yn_solid-1)+xn_solid;
        
        % Define normalized distances from marker to the down left node;
        dx_solid = zeros(size(cur_X_solid)); dy_solid = zeros(size(cur_X_solid));
%         dx_solid(X_MS_solid<gridx(nxt+1))      =   (X_MS_solid(X_MS_solid<gridx(nxt+1))-gridx(xn_solid(X_MS_solid<gridx(nxt+1)))')./xstp;
        dx_solid(X_MS_solid<gridx(nxt+1))      =   (X_MS_solid(X_MS_solid<gridx(nxt+1))-gridx(xn_solid(X_MS_solid<gridx(nxt+1)))')./(gridx(xn_solid+1)-gridx(xn_solid))';
        if nxt_aux > 0
        dx_solid(X_MS_solid>gridx(nxt+1))      =   (X_MS_solid(X_MS_solid>gridx(nxt+1))-gridx(xn_solid(X_MS_solid>gridx(nxt+1)))')./xstp_aux;
        end

%         dy_solid(Y_MS_solid<gridy(nyt+1))      =   (Y_MS_solid(Y_MS_solid<gridy(nyt+1))-gridy(yn_solid(Y_MS_solid<gridy(nyt+1)))')./ystp;
        dy_solid(Y_MS_solid<gridy(nyt+1))      =   (Y_MS_solid(Y_MS_solid<gridy(nyt+1))-gridy(yn_solid(Y_MS_solid<gridy(nyt+1)))')./(gridy(yn_solid+1)-gridy(yn_solid))';
        if nyt_aux > 0
        dy_solid(Y_MS_solid>gridy(nyt+1))      =   (Y_MS_solid(Y_MS_solid>gridy(nyt+1))-gridy(yn_solid(Y_MS_solid>gridy(nyt+1)))')./ystp_aux;
        end

        MVx_solid(:,rk)=((1.0-dx_solid).*(1.0-dy_solid).*NVx_solid(T(ielem_solid,1)) + ...
            dx_solid.*(1.0-dy_solid).*NVx_solid(T(ielem_solid,2)) + ...
            dx_solid.*dy_solid.*NVx_solid(T(ielem_solid,3)) + ...
            (1.0-dx_solid).*dy_solid.*NVx_solid(T(ielem_solid,4)));
        MVy_solid(:,rk)=((1.0-dx_solid).*(1.0-dy_solid).*NVy_solid(T(ielem_solid,1)) + ...
            dx_solid.*(1.0-dy_solid).*NVy_solid(T(ielem_solid,2)) + ...
            dx_solid.*dy_solid.*NVy_solid(T(ielem_solid,3)) + ...
            (1.0-dx_solid).*dy_solid.*NVy_solid(T(ielem_solid,4)));
        MDiv_solid(:,rk)=((1.0-dx_solid).*(1.0-dy_solid).*NDiv_solid(T(ielem_solid,1)) + ...
            dx_solid.*(1.0-dy_solid).*NDiv_solid(T(ielem_solid,2)) + ...
            dx_solid.*dy_solid.*NDiv_solid(T(ielem_solid,3)) + ...
            (1.0-dx_solid).*dy_solid.*NDiv_solid(T(ielem_solid,4)));
        
        if(rk<4)
            if (rk<3)
                cur_X_solid = X_MS_solid+dt/2*MVx_solid(:,rk);
                cur_Y_solid = Y_MS_solid+dt/2*MVy_solid(:,rk);
            else
                cur_X_solid = X_MS_solid+dt*MVx_solid(:,rk);
                cur_Y_solid = Y_MS_solid+dt*MVy_solid(:,rk);
            end
        end
        
%         xn_fluid = zeros(size(cur_X_fluid));
%         xn_fluid(X_MS_fluid<gridx(nxt+1))=double(int16(X_MS_fluid(X_MS_fluid<gridx(nxt+1))./xstp-0.5))+1;
        xn_fluid = max(cumsum(X_MS_fluid./gridx>1,2),[],2);
        if nxt_aux > 0
            xstp_aux = gridx(end)-gridx(end-1);
            xn_fluid(X_MS_fluid>gridx(nxt+1))=double(int16((X_MS_fluid(X_MS_fluid>gridx(nxt+1))-gridx(nxt+1))./xstp_aux-0.5))+1+nxt;
        end
%         yn_fluid=double(int16(Y_MS_fluid./ystp-0.5))+1;
%         yn_fluid(Y_MS_fluid<gridy(nyt+1))=double(int16(Y_MS_fluid(Y_MS_fluid<gridy(nyt+1))./ystp-0.5))+1;
        yn_fluid = max(cumsum(Y_MS_fluid./gridy>1,2),[],2);
        if nyt_aux > 0
            ystp_aux = gridy(end)-gridy(end-1);
            yn_fluid(Y_MS_fluid>gridy(nyt+1))=double(int16((Y_MS_fluid(Y_MS_fluid>gridy(nyt+1))-gridy(nyt+1))./ystp_aux-0.5))+1+nyt;
        end
        if (xn_fluid<1)
            xn_fluid(xn_fluid<1)  =   1;
        end
        if any(xn_fluid>(nxt+nxt_aux))
            xn_fluid(xn_fluid>(nxt+nxt_aux))  =   (nxt+nxt_aux);
        end
        if any(yn_fluid<1)
            yn_fluid(yn_fluid<1)  =   1;
        end
        if any(yn_fluid>(nyt+nyt_aux))
            yn_fluid(yn_fluid>(nyt+nyt_aux))  =   (nyt+nyt_aux);
        end
        % Element en el que esta ielem = (xnum-1)*(yn-1)+xn
        ielem_fluid = (nxt+nxt_aux)*(yn_fluid-1)+xn_fluid;
        
        % Define normalized distances from marker to the down left node;
        dx_fluid = zeros(size(cur_X_fluid)); dy_fluid = zeros(size(cur_X_fluid));
%         dx_fluid(X_MS_fluid<gridx(nxt+1))      =   (X_MS_fluid(X_MS_fluid<gridx(nxt+1))-gridx(xn_fluid(X_MS_fluid<gridx(nxt+1)))')./xstp;
        dx_fluid(X_MS_fluid<gridx(nxt+1))      =   (X_MS_fluid(X_MS_fluid<gridx(nxt+1))-gridx(xn_fluid(X_MS_fluid<gridx(nxt+1)))')./(gridx(xn_fluid+1)-gridx(xn_fluid))';
        if nxt_aux > 0
        dx_fluid(X_MS_fluid>gridx(nxt+1))      =   (X_MS_fluid(X_MS_fluid>gridx(nxt+1))-gridx(xn_fluid(X_MS_fluid>gridx(nxt+1)))')./xstp_aux;
        end

%         dy_fluid(Y_MS_fluid<gridy(nyt+1))      =   (Y_MS_fluid(Y_MS_fluid<gridy(nyt+1))-gridy(yn_fluid(Y_MS_fluid<gridy(nyt+1)))')./ystp;
        dy_fluid(Y_MS_fluid<gridy(nyt+1))      =   (Y_MS_fluid(Y_MS_fluid<gridy(nyt+1))-gridy(yn_fluid(Y_MS_fluid<gridy(nyt+1)))')./(gridy(yn_fluid(Y_MS_fluid<gridy(nyt+1))+1)-gridy(yn_fluid(Y_MS_fluid<gridy(nyt+1))))';
        if nyt_aux > 0
        dy_fluid(Y_MS_fluid>gridy(nyt+1))      =   (Y_MS_fluid(Y_MS_fluid>gridy(nyt+1))-gridy(yn_fluid(Y_MS_fluid>gridy(nyt+1)))')./ystp_aux;
        end

        MVx_fluid(:,rk)=((1.0-dx_fluid).*(1.0-dy_fluid).*NVx_fluid(T(ielem_fluid,1)) + ...
            dx_fluid.*(1.0-dy_fluid).*NVx_fluid(T(ielem_fluid,2)) + ...
            dx_fluid.*dy_fluid.*NVx_fluid(T(ielem_fluid,3)) + ...
            (1.0-dx_fluid).*dy_fluid.*NVx_fluid(T(ielem_fluid,4)));

        MVy_fluid(:,rk)=((1.0-dx_fluid).*(1.0-dy_fluid).*NVy_fluid(T(ielem_fluid,1)) + ...
            dx_fluid.*(1.0-dy_fluid).*NVy_fluid(T(ielem_fluid,2)) + ...
            dx_fluid.*dy_fluid.*NVy_fluid(T(ielem_fluid,3)) + ...
            (1.0-dx_fluid).*dy_fluid.*NVy_fluid(T(ielem_fluid,4)));
                
        if(rk<4)
            if (rk<3)
                cur_X_fluid = X_MS_fluid+dt/2*MVx_fluid(:,rk);
                cur_Y_fluid = Y_MS_fluid+dt/2*MVy_fluid(:,rk);
            else
                cur_X_fluid = X_MS_fluid+dt*MVx_fluid(:,rk);
                cur_Y_fluid = Y_MS_fluid+dt*MVy_fluid(:,rk);
            end
        end
        
        
    end
    
    % Recompute velocity using 4-th order Runge_Kutta
    if (markmove==4)
        MVx_solid(:,1)=(MVx_solid(:,1)+2*MVx_solid(:,2)+2*MVx_solid(:,3)+MVx_solid(:,4))/6;
        MVy_solid(:,1)=(MVy_solid(:,1)+2*MVy_solid(:,2)+2*MVy_solid(:,3)+MVy_solid(:,4))/6;
        MDiv_solid(:,1)=(MDiv_solid(:,1)+2*MDiv_solid(:,2)+2*MDiv_solid(:,3)+MDiv_solid(:,4))/6;
    end
    % Recompute velocity using 4-th order Runge_Kutta
    if (markmove==4)
        MVx_fluid(:,1)=(MVx_fluid(:,1)+2*MVx_fluid(:,2)+2*MVx_fluid(:,3)+MVx_fluid(:,4))/6;
        MVy_fluid(:,1)=(MVy_fluid(:,1)+2*MVy_fluid(:,2)+2*MVy_fluid(:,3)+MVy_fluid(:,4))/6;
    end
    
    % Displacing Marker according to its velocity
    MR_solid(:,1) = X_MS_solid + dt*MVx_solid(:,1);
%     aux_correct = MR_solid(:,1)>gridx(end) | MR_solid(:,1)<gridx(1);
%     MR_solid(aux_correct,1) = X_MS_solid(aux_correct);
    MR_solid(:,2) = Y_MS_solid + dt*MVy_solid(:,1);
%     aux_correct = MR_solid(:,2)>gridy(end) | MR_solid(:,2)<gridy(1);
%     MR_solid(aux_correct,2) = Y_MS_solid(aux_correct);
    
    % Displacing Marker according to its velocity
    MR_fluid(:,1) = X_MS_fluid + dt*MVx_fluid(:,1);
%     aux_correct = MR_fluid(:,1)>gridx(end) | MR_fluid(:,1)<gridx(1);
%     MR_fluid(aux_correct,1) = X_MS_fluid(aux_correct);
    MR_fluid(:,2) = Y_MS_fluid + dt*MVy_fluid(:,1);
%     aux_correct = MR_fluid(:,2)>gridy(end) | MR_fluid(:,2)<gridy(1);
%     MR_fluid(aux_correct,2) = Y_MS_fluid(aux_correct);
    
    dt_poros = dt;
    
    % Solve n
    a = repmat(MDiv_solid(:,1),1,nTP);
    c = MR_solid(:,end-nTP+1:end);
    
    MR_solid(:,end-nTP+1:end) = c.*exp(-a.*dt_poros);
    
end



                