function test1

P = 20; % number of channels
Q = 3;  % model order

a_true = randn(P, Q);
b_true = randn(Q, Q, Q)+sqrt(-1)*randn(Q, Q, Q);
b_true = b_true./norm(b_true(:));

B = tprod(a_true, [3 -1], tprod(a_true, [2 -1], tprod(a_true, [1 -1], b_true, [-1 2 3]), [1 -1 3]), [1 2 -1]);
disp(size(B))


a_init = randn(P, Q);
b_init = randn(Q, Q, Q);
x0 = [vec(a_init); vec(b_init)];
% x0 = [vec(a_true); vec(b_true)];

opts = optimset('Algorithm','levenberg-marquardt','Display','on', 'Jacobian', 'on', 'DerivativeCheck', 'off');

tic
[x resnorm] = lsqnonlin(@(x) myfun(x), x0, [], [], opts);
toc

a_est = real(reshape(x(1:P*Q), P, Q));
b_est = reshape(x((P*Q+1):end), [Q Q Q]);

norm(vec(B-tprod(a_est, [3 -1], tprod(a_est, [2 -1], tprod(a_est, [1 -1], b_est, [-1 2 3]), [1 -1 3]), [1 2 -1])))/norm(vec(B))

subspace(a_true, a_est)

keyboard

function [F J] = myfun(x)

a = real(reshape(x(1:P*Q), P, Q));
b = reshape(x((P*Q+1):end), [Q Q Q]);

F = (vec(B-tprod(a, [3 -1], tprod(a, [2 -1], tprod(a, [1 -1], b, [-1 2 3]), [1 -1 3]), [1 2 -1])));

J = zeros(P, P, P, P, Q);

for qq = 1:Q
  aba1(:, :, qq) = a*squeeze(b(qq, :, :))*a';
  aba2(:, :, qq) = a*squeeze(b(:, qq, :))*a';
  aba3(:, :, qq) = a*squeeze(b(:, :, qq))*a';
end

for pp = 1:P
  J(pp, :, :, pp, :) = squeeze(J(pp, :, :, pp, :)) - aba1;
  J(:, pp, :, pp, :) = squeeze(J(:, pp, :, pp, :)) - aba2;
  J(:, :, pp, pp, :) = squeeze(J(:, :, pp, pp, :)) - aba3;
end

J = [reshape(J, P^3, P*Q) -kron(kron(a, a), a)];

end

end

