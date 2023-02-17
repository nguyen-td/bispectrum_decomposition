function test1

P = 60;
Q = 3;

snr = 0.75;
nstart = 1;
i = sqrt(-1);

a_true = randn(P, Q);
b_true = randn(Q, Q, Q)+i*randn(Q, Q, Q);

B = tprod(a_true, [3 -1], tprod(a_true, [2 -1], tprod(a_true, [1 -1], b_true, [-1 2 3]), [1 -1 3]), [1 2 -1]);
no = randn(size(B)) + i*randn(size(B));
no = norm(vec(B))*no/norm(vec(no));
B = snr*B + (1-snr)*no;

for istart = 1:nstart
  a_init(:, :, istart) = randn(P, Q);
  b_init(:, :, :, istart) = randn(Q, Q, Q) + i*randn(Q, Q, Q);
end

% levenberg marquardt
opts = optimset('Algorithm','levenberg-marquardt','Display','off', 'Jacobian', 'on', 'FinDiffType', 'central', 'DerivativeCheck', 'off');

tic
for istart = 1:nstart
  disp(istart)
  [x_lm(:, istart) resnorm_lm(istart)] = lsqnonlin(@(x) myfun_lm(x), [vec(a_init(:, :, istart)); vec(b_init(:, :, :, istart))], [], [], opts);
end
toc

[f_lm in] = min(resnorm_lm);
x = x_lm(:, in);
a_lm = real(reshape(x(1:P*Q), P, Q));
b_lm = reshape(x((P*Q+1):end), [Q Q Q]);

f_lm
angle_lm = subspace(a_true, a_lm)


% l-BFGS
opts = optimset('Display','off', 'GradObj', 'on', 'FinDiffType', 'central', 'DerivativeCheck', 'off');

tic
for istart = 1:nstart
  disp(istart)
  [x_lbfgs(:, istart) resnorm_lbfgs(istart)] = fminunc(@(x) myfun_lbfgs(x), [vec(a_init(:, :, istart)); vec(real(b_init(:, :, :, istart))); vec(imag(b_init(:, :, :, istart)))], opts);
end
toc

[f_lbfgs in] = min(resnorm_lbfgs);
x = x_lbfgs(:, in);
a_lbfgs = reshape(x(1:P*Q), P, Q);
b_lbfgs = reshape(x((P*Q+1):end), [Q Q Q 2]);
b_lbfgs = b_lbfgs(:, :, :, 1) + i*b_lbfgs(:, :, :, 2);

f_lbfgs
angle_lbfgs = subspace(a_true, a_lbfgs)


keyboard

% computes unsquared errors and Jacobian for levenberg marquard
function [f J] = myfun_lm(x)

  a = real(reshape(x(1:P*Q), P, Q));
  b = reshape(x((P*Q+1):end), [Q Q Q]);

  f = (vec(B-tprod(a, [3 -1], tprod(a, [2 -1], tprod(a, [1 -1], b, [-1 2 3]), [1 -1 3]), [1 2 -1])));

  J = zeros(P, P, P, P, Q);

  for qq = 1:Q
    aba1(:, :, qq) = a*squeeze(b(qq, :, :))*a';
    aba2(:, :, qq) = a*squeeze(b(:, qq, :))*a';
    aba3(:, :, qq) = a*squeeze(b(:, :, qq))*a';
  end

  for pp = 1:P
    J(pp, :, :, pp, :) = -aba1;
    J(:, pp, :, pp, :) = squeeze(J(:, pp, :, pp, :)) - aba2;
    J(:, :, pp, pp, :) = squeeze(J(:, :, pp, pp, :)) - aba3;
  end
  
  J = [reshape(J, P^3, P*Q) -kron(kron(a, a), a)];
  
end

% computes summed squared errors and gradient for lBFGS algorithm
% memory efficient implementation
function [f grad] = myfun_lbfgs(x)

  a = reshape(x(1:P*Q), P, Q);
  b = reshape(x((P*Q+1):end), [Q Q Q 2]);
  b = b(:, :, :, 1) + i*b(:, :, :, 2);

  f = vec(B-tprod(a, [3 -1], tprod(a, [2 -1], tprod(a, [1 -1], b, [-1 2 3]), [1 -1 3]), [1 2 -1]));

  J = zeros(P, P, P, P, Q);
  for qq = 1:Q
    aba1(:, :, qq) = a*squeeze(b(qq, :, :))*a';
    aba2(:, :, qq) = a*squeeze(b(:, qq, :))*a';
    aba3(:, :, qq) = a*squeeze(b(:, :, qq))*a';
  end
  
  grada = zeros(P, Q);
  for pp = 1:P
    for qq = 1:Q
      J = zeros(P, P, P);
      J(pp, :, :) = -aba1(:, :, qq);
      J(:, pp, :) = squeeze(J(:, pp, :)) - aba2(:, :, qq);
      J(:, :, pp) = squeeze(J(:, :, pp)) - aba3(:, :, qq);
      grada(pp, qq) = 2*real(f)'*real(vec(J)) + 2*imag(f)'*imag(vec(J));
    end
  end
  grada = vec(grada);
  
  gradb = zeros(Q, Q, Q);
  for q1 = 1:Q
    for q2 = 1:Q
      for q3 = 1:Q
        gradb(q1, q2, q3) = 2*f'*kron(kron(a(:, q3), a(:, q2)), a(:, q1));
      end
    end
  end
  gradb = vec(gradb);
  
  grad = [grada; -real(gradb); imag(gradb)];

  f = sum(real(f).^2 + imag(f).^2);    
end


% computes summed squared errors and gradient for lBFGS algorithm
% creates full Jacobian as intermediate step - memory inefficient
function [f grad] = myfun_lbfgs2(x)

  a = reshape(x(1:P*Q), P, Q);
  b = reshape(x((P*Q+1):end), [Q Q Q 2]);
  b = b(:, :, :, 1) + i*b(:, :, :, 2);

  f = vec(B-tprod(a, [3 -1], tprod(a, [2 -1], tprod(a, [1 -1], b, [-1 2 3]), [1 -1 3]), [1 2 -1]));

  J = zeros(P, P, P, P, Q);

  for qq = 1:Q
    aba1(:, :, qq) = a*squeeze(b(qq, :, :))*a';
    aba2(:, :, qq) = a*squeeze(b(:, qq, :))*a';
    aba3(:, :, qq) = a*squeeze(b(:, :, qq))*a';
  end

  for pp = 1:P
    J(pp, :, :, pp, :) = -aba1;
    J(:, pp, :, pp, :) = squeeze(J(:, pp, :, pp, :)) - aba2;
    J(:, :, pp, pp, :) = squeeze(J(:, :, pp, pp, :)) - aba3;
  end
  
  Ja = reshape(J, P^3, P*Q);
  
  grada = 2*real(f)'*real(Ja) + 2*imag(f)'*imag(Ja);
  gradb = 2*f'*kron(kron(a, a), a);
  grad = [grada -real(gradb) imag(gradb)]';

  f = sum(real(f).^2 + imag(f).^2);    
end

function v = vec(x)
  v = reshape(x, [], 1);  
end

end
