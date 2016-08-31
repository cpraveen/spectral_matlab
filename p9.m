% p9.m - polynomial interpolation in equispaced and Chebyshev pts

  N = 16;
  xx = -1.01:.005:1.01; clf
  for i = 1:2
    if i==1, s = 'equispaced points'; x = -1 + 2*(0:N)/N; end
    if i==2, s = 'Chebyshev points';  x = cos(pi*(0:N)/N); end
    subplot(2,2,i)
    u = 1./(1+16*x.^2);
    uu = 1./(1+16*xx.^2);
    p = polyfit(x,u,N);              % interpolation
    pp = polyval(p,xx);              % evaluation of interpolant
    plot(x,u,'.','markersize',13)
    line(xx,pp)
    axis([-1.1 1.1 -1 1.5]), title(s)
    error = norm(uu-pp,inf);
    text(-.5,-.5,['max error = ' num2str(error)])
  end
