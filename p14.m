% p14.m - solve nonlinear BVP u_xx = exp(u), u(-1)=u(1)=0
%         (compare p13.m)

  N = 16;
  [D,x] = cheb(N); D2 = D^2; D2 = D2(2:N,2:N);
  u = zeros(N-1,1);
  change = 1; it = 0;
  while change > 1e-15                   % fixed-point iteration
    unew = D2\exp(u); 
    change = norm(unew-u,inf);
    u = unew; it = it+1;
  end
  u = [0;u;0];
  clf, subplot('position',[.1 .4 .8 .5])
  plot(x,u,'.','markersize',16)
  xx = -1:.01:1;
  uu = polyval(polyfit(x,u,N),xx);
  line(xx,uu), grid on
  title(sprintf('no. steps = %d      u(0) =%18.14f',it,u(N/2+1)))
