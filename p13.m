% p13.m - solve linear BVP u_xx = exp(4x), u(-1)=u(1)=0

  N = 16;
  [D,x] = cheb(N);        
  D2 = D^2;              
  D2 = D2(2:N,2:N);                   % boundary conditions
  f = exp(4*x(2:N));           
  u = D2\f;                           % Poisson eq. solved here
  u = [0;u;0];              
  clf, subplot('position',[.1 .4 .8 .5])
  plot(x,u,'.','markersize',16)
  xx = -1:.01:1;
  uu = polyval(polyfit(x,u,N),xx);    % interpolate grid data
  line(xx,uu)
  grid on
  exact = ( exp(4*xx) - sinh(4)*xx - cosh(4) )/16; 
  title(['max err = ' num2str(norm(uu-exact,inf))],'fontsize',12)
