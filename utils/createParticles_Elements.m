function [MS_C_solid,MS_C_fluid] = createParticles_Elements(mxelem,myelem,nxv,nyv,xsize,ysize,nTE,nTP)
% Create the particles and storage room

mxnum = 1;                   % total number of markers in horizontal direction
mynum_s = 2*nyv*myelem;      % total number of markers in vertical direction
mynum_f = 2*nyv*myelem;      % total number of markers in vertical direction
mxstep = xsize/mxnum;        % step between markers in horizontal direction
mystep_s = ysize/mynum_s;    % step between markers in vertical direction
mystep_f = ysize/mynum_f;    % step between markers in vertical direction

% Particle Structure in a matrix
%                                # Properties
%                             ------------------->
%                   X - Y - phi - C_TE x nTE - n_TP x nTP
%            |  1
%  # Markers |  2
%            |  3
%            v  4

MS_aux = sparse(mynum_s*mxnum,2);
mm1 = 0;         % marker counter

stringx = mxstep/2:mxstep:mxnum*mxstep-mxstep/2;
stringx = 0;
MX_s = repmat(stringx,mynum_s,1);
MX_f = repmat(stringx,mynum_f,1);
MX_rand_s = (rand(size(MX_s))-0.5)*mxstep;
MX_rand_f = (rand(size(MX_f))-0.5)*mxstep;
stringy_s = mystep_s/2:mystep_s:mynum_s*mystep_s-mystep_s/2;
stringy_f = mystep_f/2:mystep_f:mynum_f*mystep_f-mystep_f/2;
MY_s = repmat(stringy_s',1,mxnum);
MY_f = repmat(stringy_f',1,mxnum);
MY_rand_s = (rand(size(MY_s))-0.5)*mystep_s;
MY_rand_f = (rand(size(MY_f))-0.5)*mystep_f;
MS_aux_s = [MX_s(:)+0*MX_rand_s(:) MY_s(:)+0*MY_rand_s(:)];
MS_aux_f = [MX_f(:)+0*MX_rand_f(:) MY_f(:)+0*MY_rand_f(:)];
% MS_aux = [MX(:)+0*MX_rand(:) MY(:)+0*MY_rand(:)];

MS_C_s = zeros(mynum_s*mxnum,3+nTE+nTP);
MS_C_f = zeros(mynum_f*mxnum,3+nTE+nTP);
MS_C_s(:,[1 2]) = MS_aux_s;
MS_C_f(:,[1 2]) = MS_aux_f;
MS_C_solid = MS_C_s;
MS_C_fluid = MS_C_f;

