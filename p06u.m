% p6u.m - variable coefficient wave equation - UNSTABLE VARIANT

% Grid, variable coefficient, and initial data:
  N = 128; h = 2*pi/N; x = h*(1:N);
  c = .2 + sin(x-1).^2;
  t = 0; dt = 1.9/N;
  v = exp(-100*(x-1).^2);
  vold = exp(-100*(x-.2*dt-1).^2);

% Time-stepping by leap frog formula:
  tmax = 8; tplot = .15; clf, set(gcf,'renderer','zbuffer')
  plotgap = round(tplot/dt);
  nplots = round(tmax/tplot);
  data = [v; zeros(nplots,N)]; tdata = t;
  for i = 1:nplots
    for n = 1:plotgap
      t = t+dt;
      v_hat = fft(v);
      w_hat = 1i*[0:N/2-1 0 -N/2+1:-1] .* v_hat;
      w = real(ifft(w_hat)); 
      vnew = vold - 2*dt*c.*w;        % leap frog formula
      vold = v; v = vnew;
    end
    data(i+1,:) = v; tdata = [tdata; t];
    if max(abs(v)>2.5), data = data(1:i+1,:); break, end
  end

% Plot results:
  waterfall(x,tdata,data)
  axis([0 2*pi 0 tmax -3. 3.]), view(10,70), grid off
  colormap(1e-6*[1 1 1]); xlabel x, ylabel t, zlabel u
