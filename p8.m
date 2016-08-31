% p8.m - eigenvalues of harmonic oscillator -u"+x^2 u on R

  format long, format compact
  L = 8;                             % domain is [-L L], periodic
  for N = 6:6:36
    h = 2*pi/N; x = h*(1:N); x = L*(x-pi)/pi;
    column = [-pi^2/(3*h^2)-1/6 ...
       -.5*(-1).^(1:N-1)./sin(h*(1:N-1)/2).^2];
    D2 = (pi/L)^2*toeplitz(column);  % 2nd-order differentiation
    eigenvalues = sort(eig(-D2 + diag(x.^2)));
    N, eigenvalues(1:4)
  end
