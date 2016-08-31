% p36.m - Laplace eq. on [-1,1]x[-1,1] with nonzero BCs

% Set up grid and 2D Laplacian, boundary points included:
  N = 24; [D,x] = cheb(N); y = x;
  [xx,yy] = meshgrid(x,y); xx = xx(:); yy = yy(:);        
  D2 = D^2; I = eye(N+1); L = kron(I,D2) + kron(D2,I);     

% Impose boundary conditions by replacing appropriate rows of L:
  b = find(abs(xx)==1 | abs(yy)==1);            % boundary pts
  L(b,:) = zeros(4*N,(N+1)^2); L(b,b) = eye(4*N);
  rhs = zeros((N+1)^2,1);
  rhs(b) = (yy(b)==1).*(xx(b)<0).*sin(pi*xx(b)).^4 + ...
           .2*(xx(b)==1).*sin(3*pi*yy(b));

% Solve Laplace equation, reshape to 2D, and plot:
  u = L\rhs; uu = reshape(u,N+1,N+1);
  [xx,yy] = meshgrid(x,y);
  [xxx,yyy] = meshgrid(-1:.04:1,-1:.04:1);
  uuu = interp2(xx,yy,uu,xxx,yyy,'cubic');
  clf, subplot('position',[.1 .4 .8 .5])
  mesh(xxx,yyy,uuu), colormap(1e-6*[1 1 1]);
  axis([-1 1 -1 1 -.2 1]), view(-20,45)
  text(0,.8,.4,sprintf('u(0,0) = %12.10f',uu(N/2+1,N/2+1)))
