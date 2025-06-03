function x = solveWithLagrangeMultipliers( A, b, C, method )
% x = solveWithLagrangeMultipliers( A, b, C, method )
%
% Solve the system Ax=b with the boundary conditions decribed in C 
% using Lagrange multipliers
% 
% INPUT
%     A       matrix
%     b       vector
%     C       boundary conditions as used in useSBC_0.m. size(C)=(n,2)
%             C(:,1) = positions in the matrix of the prescribed values
%             C(:,2) = prescribed values
%     method  1=bigc, 0= \
%
% OUTPUT
%     x     solution
%
if nargin < 4
   method = 0;
end
if ~isempty( C )
   nunk = size( A, 1 );
   nDir = size( C, 1 );
   Acc = zeros( nDir, nunk );
   Acc(:,C(:,1)) = eye( nDir );
   bcc = C(:,2);
   Z0 = zeros( nDir, nDir );
   Atot = [A    Acc'; ...
           Acc  Z0  ];
   btot = [b; bcc];
   if method
      [x,flag,relres,iter,resvec] = bicg( Atot, btot );
   else
      x = Atot\btot;
   end
   x = x(1:nunk);
else
   if method
      [x,flag,relres,iter,resvec] = gmres( A, b );
   else
      x = A\b;
   end
end



