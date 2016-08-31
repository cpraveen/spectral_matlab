% p7.m - accuracy of periodic spectral differentiation

% Compute derivatives for various values of N:
  Nmax = 50; E = zeros(3,Nmax/2-2);
  for N = 6:2:Nmax;
    h = 2*pi/N; x = h*(1:N)';
    column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)]';
    D = toeplitz(column,column([1 N:-1:2]));
    v = abs(sin(x)).^3;                     % 3rd deriv in BV
    vprime = 3*sin(x).*cos(x).*abs(sin(x));
    E(1,N/2-2) = norm(D*v-vprime,inf);
    v = exp(-sin(x/2).^(-2));               % C-infinity
    vprime = .5*v.*sin(x)./sin(x/2).^4;
    E(2,N/2-2) = norm(D*v-vprime,inf);
    v = 1./(1+sin(x/2).^2);                 % analytic in a strip
    vprime = -sin(x/2).*cos(x/2).*v.^2;
    E(3,N/2-2) = norm(D*v-vprime,inf);
    v = sin(10*x); vprime = 10*cos(10*x);   % band-limited
    E(4,N/2-2) = norm(D*v-vprime,inf);
  end

% Plot results:
  titles = {'|sin(x)|^3','exp(-sin^{-2}(x/2))',...
       '1/(1+sin^2(x/2))','sin(10x)'}; clf
  for iplot = 1:4 
    subplot(2,2,iplot)
    semilogy(6:2:Nmax,E(iplot,:),'.','markersize',12)
    line(6:2:Nmax,E(iplot,:))
    axis([0 Nmax 1e-16 1e3]), grid on
    set(gca,'xtick',0:10:Nmax,'ytick',(10).^(-15:5:0))
    xlabel N, ylabel error, title(titles(iplot))
  end
