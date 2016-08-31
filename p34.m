% p34.m - Allen-Cahn eq. u_t = eps*u_xx+u-u^3, u(-1)=-1, u(1)=1
%         (compare p6.m and p32.m)

% Differentiation matrix and initial data:
  N = 20; [D,x] = cheb(N); D2 = D^2;     % use full-size matrix
  D2([1 N+1],:) = zeros(2,N+1);          % for convenience
  eps = 0.01; dt = min([.01,50*N^(-4)/eps]);
  t = 0; v = .53*x + .47*sin(-1.5*pi*x);

% Solve PDE by Euler formula and plot results:
  tmax = 100; tplot = 2; nplots = round(tmax/tplot);
  plotgap = round(tplot/dt); dt = tplot/plotgap;
  xx = -1:.025:1; vv = polyval(polyfit(x,v,N),xx);
  plotdata = [vv; zeros(nplots,length(xx))]; tdata = t;
  for i = 1:nplots
    for n = 1:plotgap
      t = t+dt; v = v + dt*(eps*D2*(v-x) + v - v.^3);    % Euler
    end
    vv = polyval(polyfit(x,v,N),xx);
    plotdata(i+1,:) = vv; tdata = [tdata; t];
  end
  clf, subplot('position',[.1 .4 .8 .5])
  mesh(xx,tdata,plotdata), grid on, axis([-1 1 0 tmax -1 1]),
  view(-60,55), colormap(1e-6*[1 1 1]); xlabel x, ylabel t, zlabel u
