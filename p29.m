% p29.m - solve Poisson equation on the unit disk
%         (compare p16.m and p28.m)

% Laplacian in polar coordinates:
  N = 31; [D,r] = cheb(N); N2 = (N-1)/2; D2 = D^2;
  D1 = D2(2:N2+1,2:N2+1); D2 = D2(2:N2+1,N:-1:N2+2);
  E1 =  D(2:N2+1,2:N2+1); E2 =  D(2:N2+1,N:-1:N2+2);
  M = 40; dt = 2*pi/M; t = dt*(1:M)'; M2 = M/2;
  D2t = toeplitz([-pi^2/(3*dt^2)-1/6 ...
                   .5*(-1).^(2:M)./sin(dt*(1:M-1)/2).^2]);
  R = diag(1./r(2:N2+1)); Z = zeros(M2); I = eye(M2);
  L = kron(D1+R*E1,eye(M))+kron(D2+R*E2,[Z I;I Z])+kron(R^2,D2t);

% Right-hand side and solution for u:
  [rr,tt] = meshgrid(r(2:N2+1),t); rr = rr(:); tt = tt(:);
  f = -rr.^2.*sin(tt/2).^4 + sin(6*tt).*cos(tt/2).^2; u = L\f;

% Reshape results onto 2D grid and plot them:
  u = reshape(u,M,N2); u = [zeros(M+1,1) u([M 1:M],:)];
  [rr,tt] = meshgrid(r(1:N2+1),t([M 1:M]));
  [xx,yy] = pol2cart(tt,rr);
  clf, subplot('position',[.1 .4 .8 .5])
  mesh(xx,yy,u), view(20,40), colormap(1e-6*[1 1 1]);
  axis([-1 1 -1 1 -.01 .05]), xlabel x, ylabel y, zlabel u 
