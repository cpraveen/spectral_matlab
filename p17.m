% p17.m - Helmholtz eq. u_xx + u_yy + (k^2)u = f
%         on [-1,1]x[-1,1]    (compare p16.m)

% Set up spectral grid and tensor product Helmholtz operator:
  N = 24; [D,x] = cheb(N); y = x;
  [xx,yy] = meshgrid(x(2:N),y(2:N));
  xx = xx(:); yy = yy(:);
  f = exp(-10*((yy-1).^2+(xx-.5).^2));
  D2 = D^2; D2 = D2(2:N,2:N); I = eye(N-1);
  k = 9;
  L = kron(I,D2) + kron(D2,I) + k^2*eye((N-1)^2);

% Solve for u, reshape to 2D grid, and plot:
  u = L\f;
  uu = zeros(N+1,N+1); uu(2:N,2:N) = reshape(u,N-1,N-1);
  [xx,yy] = meshgrid(x,y);
  [xxx,yyy] = meshgrid(-1:.0333:1,-1:.0333:1);
  uuu = interp2(xx,yy,uu,xxx,yyy,'cubic');
  figure(1), clf, mesh(xxx,yyy,uuu), colormap(1e-6*[1 1 1]);
  xlabel x, ylabel y, zlabel u
  text(.2,1,.022,sprintf('u(0,0) = %13.11f',uu(N/2+1,N/2+1)))
  figure(2), clf, contour(xxx,yyy,uuu)
  colormap(1e-6*[1 1 1]); axis square
