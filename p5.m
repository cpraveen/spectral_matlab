% p5.m - repetition of p4.m via FFT
%        For complex v, delete "real" commands.

% Differentiation of a hat function:
  N = 24; h = 2*pi/N; x = h*(1:N)';
  v = max(0,1-abs(x-pi)/2); v_hat = fft(v);
  w_hat = 1i*[0:N/2-1 0 -N/2+1:-1]' .* v_hat;
  w = real(ifft(w_hat)); clf
  subplot(3,2,1), plot(x,v,'.-','markersize',13)
  axis([0 2*pi -.5 1.5]), grid on, title('function')
  subplot(3,2,2), plot(x,w,'.-','markersize',13)
  axis([0 2*pi -1 1]), grid on, title('spectral derivative')

% Differentiation of exp(sin(x)):
  v = exp(sin(x)); vprime = cos(x).*v;
  v_hat = fft(v);
  w_hat = 1i*[0:N/2-1 0 -N/2+1:-1]' .* v_hat;
  w = real(ifft(w_hat));
  subplot(3,2,3), plot(x,v,'.-','markersize',13)
  axis([0 2*pi 0 3]), grid on
  subplot(3,2,4), plot(x,w,'.-','markersize',13)
  axis([0 2*pi -2 2]), grid on
  error = norm(w-vprime,inf);
  text(2.2,1.4,['max error = ' num2str(error)])
