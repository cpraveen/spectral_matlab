
  x = x(:); M = length(x); a = zeros(M,1); D = zeros(M,M);
  for k = 1:M; a(k) = prod(x(k)-x([1:k-1 k+1:M])); end
  for j = 1:M; t = x(:)-x(j); t(j) = 1; D(:,j) = 1./t; end
  A = diag(a); D = A*D/A;
  for j = 1:M
    sum = 0;
    for k = [1:j-1 j+1:M];
      sum = sum + 1/(x(j)-x(k));
    end
    D(j,j) = sum;
  end
