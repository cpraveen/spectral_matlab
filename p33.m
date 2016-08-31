% p33.m - solve linear BVP u_xx = exp(4x), u'(-1)=u(1)=0

  N = 16;
  [D,x] = cheb(N); D2 = D^2;              
  D2(N+1,:) = D(N+1,:);            % Neumann condition at x = -1
  D2 = D2(2:N+1,2:N+1);                
  f = exp(4*x(2:N));           
  u = D2\[f;0];                       
  u = [0;u];              
  clf, subplot('position',[.1 .4 .8 .5])
  plot(x,u,'.','markersize',16)
  axis([-1 1 -4 0])
  xx = -1:.01:1;
  uu = polyval(polyfit(x,u,N),xx);   
  line(xx,uu)
  grid on
  exact = (exp(4*xx) - 4*exp(-4)*(xx-1) - exp(4))/16; 
  title(['max err = ' num2str(norm(uu-exact,inf))],'fontsize',12)
