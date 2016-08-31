% p20u.m - 2nd-order wave eq. in 2D via FFT (compare p19.m) UNSTABLE

% Grid and initial data:
  N = 24; x = cos(pi*(0:N)/N); y = x'; dt = 6.6/N^2;
  [xx,yy] = meshgrid(x,y);
  plotgap = round((1/3)/dt);
  vv = exp(-40*((xx-.4).^2 + yy.^2));
  vvold = vv; vvnew = 0*vv;

% Time-stepping by leap frog formula:
  [ay,ax] = meshgrid([.56 .06],[.1 .55]);
  for n = 0:3*plotgap
    t = n*dt;
    if rem(n+.5,plotgap)<1       % Plots at multiples of t=1/3
      i = n/plotgap+1; subplot('position',[ax(i) ay(i) .36 .36])
      [xxx,yyy] = meshgrid(-1:.05:1,-1:.05:1);
      vvv = interp2(xx,yy,vv,xxx,yyy,'cubic');
      mesh(xxx,yyy,vvv), axis([-1 1 -1 1 -0.15 1])
      colormap(1e-6*[1 1 1]); title(['t = ' num2str(t)]), drawnow
    end
    uxx = zeros(N+1,N+1); uyy = zeros(N+1,N+1);
    ii = 2:N;
    for i = ii                 % 2nd derivs wrt x in each row
      v = vv(i,:); V = [v fliplr(v(ii))]; U = real(fft(V));
      W1 = real(ifft(1i*[0:N-1 0 1-N:-1].*U)); % diff wrt theta
      W2 = real(ifft(-[0:N 1-N:-1].^2.*U));    % diff^2 wrt theta
      uxx(i,ii) = W2(ii)./(1-x(ii).^2) - x(ii).* ... 
                     W1(ii)./(1-x(ii).^2).^(3/2);
    end
    for j = ii                 % 2nd derivs wrt y in each column
      v = vv(:,j); V = [v; flipud(v(ii))]; U = real(fft(V));
      W1 = real(ifft(1i*[0:N-1 0 1-N:-1]'.*U));% diff wrt theta   
      W2 = real(ifft(-[0:N 1-N:-1]'.^2.*U));   % diff^2 wrt theta
      uyy(ii,j) = W2(ii)./(1-y(ii).^2) - y(ii).* ...
                     W1(ii)./(1-y(ii).^2).^(3/2);
    end
    vvnew = 2*vv - vvold + dt^2*(uxx+uyy); 
    vvold = vv; vv = vvnew;
    if max(max(abs(vv)))>2, break, end
  end

close
      [xxx,yyy] = meshgrid(-1:.05:1,-1:.05:1);
      vvv = interp2(xx,yy,vv,xxx,yyy,'cubic');
      mesh(xxx,yyy,vvv), colormap(1e-6*[1 1 1]);
      axis([-1 1 -1 1 -1 1]), title(['t = ' num2str(t)]), drawnow

