% p39.m - eigenmodes of biharmonic on a square with clamped BCs
%         (compare p38.m)

% Construct spectral approximation to biharmonic operator:
  N = 17; [D,x] = cheb(N); D2 = D^2; D2 = D2(2:N,2:N);
  S = diag([0; 1 ./(1-x(2:N).^2); 0]);
  D4 = (diag(1-x.^2)*D^4 - 8*diag(x)*D^3 - 12*D^2)*S;
  D4 = D4(2:N,2:N); I = eye(N-1);
  L = kron(I,D4) + kron(D4,I) + 2*kron(D2,I)*kron(I,D2);

% Find and plot 25 eigenmodes:
  [V,Lam] = eig(-L); Lam = -real(diag(Lam));
  [Lam,ii] = sort(Lam); ii = ii(1:25); V = real(V(:,ii));
  Lam = sqrt(Lam/Lam(1));
  [xx,yy] = meshgrid(x,x); 
  [xxx,yyy] = meshgrid(-1:.01:1,-1:.01:1);
  [ay,ax] = meshgrid(.8:-.2:0,0:.16:.64);
  sq = [1+1i -1+1i -1-1i 1-1i 1+1i]; clf
  for i = 1:25
    uu = zeros(N+1,N+1); uu(2:N,2:N) = reshape(V(:,i),N-1,N-1);
    subplot('position',[ax(i) ay(i) .16 .2]), plot(sq)
    uuu = interp2(xx,yy,uu,xxx,yyy,'cubic');
    hold on, contour(xxx,yyy,uuu,[0 0]), axis square
    axis (1.25*[-1 1 -1 1]), axis off, colormap(1e-6*[1 1 1]);
    text(-.3,1.15,num2str(Lam(i)),'fontsize',7)
  end
