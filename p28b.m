% p28b.m - eigenmodes of Laplacian on the disk

% r coordinate, ranging from -1 to 1 (N must be odd)
  N = 25; [D,r] = cheb(N); N2 = (N-1)/2; D2 = D^2;
  D1 = D2(2:N2+1,2:N2+1); D2 = D2(2:N2+1,N:-1:N2+2);
  E1 =  D(2:N2+1,2:N2+1); E2 =  D(2:N2+1,N:-1:N2+2);

% t = theta coordinate, ranging from 0 to 2*pi (M must be even):
  M = 20; dt = 2*pi/M; t = dt*(1:M)'; M2 = M/2;
  D2t = toeplitz([-pi^2/(3*dt^2)-1/6 ...
                    .5*(-1).^(2:M)./sin(dt*(1:M-1)/2).^2]);

% Laplacian in polar coordinates:
  R = diag(1./r(2:N2+1)); Z = zeros(M2); I = eye(M2);
  L = kron(D1+R*E1,eye(M))+ kron(D2+R*E2,[Z I;I Z])+ kron(R^2,D2t);

% Compute four eigenmodes:
  index = 1:25;
  [V,Lam] = eig(-L); Lam = diag(Lam);
  [Lam,ii] = sort(real(Lam)); ii = ii(index); V = real(V(:,ii));
  Lam = sqrt(Lam(index)/Lam(1));

% Plot eigenmodes with nodal lines underneath:
  [rr,tt] = meshgrid(r(1:N2+1),[0;t]); [xx,yy] = pol2cart(tt,rr);
  z = exp(sqrt(-1)*pi*(-100:100)/100);
  [ay,ax] = meshgrid(.8:-.2:0,0:.16:.64); clf
  for i = 1:25
    u = reshape(real(V(:,i)),M,N2);
    u = [zeros(M+1,1) u([M 1:M],:)];
    u = u/norm(u(:),inf);
    subplot('position',[ax(i) ay(i) .16 .2]), plot(z)
    axis(1.25*[-1 1 -1 1 -1 1]), axis off, hold on
    view(0,90), colormap(1e-6*[1 1 1]); axis square
    contour3(xx,yy,u-1,[-1 -1]), plot3(real(z),imag(z),-abs(z))
    text(-.3,1.15,num2str(Lam(i)),'fontsize',7)
  end
