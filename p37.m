% p37.m - 2D "wave tank" with Neumann BCs for |y|=1

% x variable in [-A,A], Fourier:
  A = 3; Nx = 50; dx = 2*A/Nx; x = -A+dx*(1:Nx)';
  D2x = (pi/A)^2*toeplitz([-1/(3*(dx/A)^2)-1/6 ...
           .5*(-1).^(2:Nx)./sin((pi*dx/A)*(1:Nx-1)/2).^2]);

% y variable in [-1,1], Chebyshev:
  Ny = 15; [Dy,y] = cheb(Ny); D2y = Dy^2;
  BC = -Dy([1 Ny+1],[1 Ny+1])\Dy([1 Ny+1],2:Ny);

% Grid and initial data:
  [xx,yy] = meshgrid(x,y);
  vv = exp(-8*((xx+1.5).^2+yy.^2));
  dt = 5/(Nx+Ny^2); 
  vvold = exp(-8*((xx+dt+1.5).^2+yy.^2));

% Time-stepping by leap frog formula:
  clf, plotgap = round(2/dt); dt = 2/plotgap;
  for n = 0:2*plotgap
    t = n*dt;
    if rem(n+.5,plotgap)<1
      subplot(3,1,n/plotgap+1), mesh(xx,yy,vv), view(-10,60)
      axis([-A A -1 1 -0.15 1]), colormap(1e-6*[1 1 1]);
      text(-2.5,1,.5,['t = ' num2str(t)],'fontsize',18),
      set(gca,'ztick',[]), grid off, drawnow
    end
    vvnew = 2*vv - vvold + dt^2*(vv*D2x +D2y*vv);
    vvold = vv; vv = vvnew;
    vv([1 Ny+1],:) = BC*vv(2:Ny,:);       % Neumann BCs for |y|=1
  end
