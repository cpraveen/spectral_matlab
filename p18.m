% p18.m - Chebyshev differentiation via FFT (compare p11.m)

  xx = -1:.01:1; ff = exp(xx).*sin(5*xx); clf
  for N = [10 20]
    x = cos(pi*(0:N)'/N); f = exp(x).*sin(5*x);
      subplot('position',[.15 .66-.4*(N==20) .31 .28])
      plot(x,f,'.','markersize',14), grid on
      line(xx,ff)
      title(['f(x), N=' int2str(N)])
    error = chebfft(f) - exp(x).*(sin(5*x)+5*cos(5*x));
      subplot('position',[.55 .66-.4*(N==20) .31 .28])
      plot(x,error,'.','markersize',14), grid on
      line(x,error)
      title(['error in f''(x), N=' int2str(N)])
  end

  
