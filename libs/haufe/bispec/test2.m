function test1

P = 10;
Q = 3;

a_true = randn(P, Q);
b_true = randn(Q, Q, Q)+sqrt(-1)*randn(Q, Q, Q);
b_true = b_true./norm(vec(b_true));

B = tprod(a_true, [3 -1], tprod(a_true, [2 -1], tprod(a_true, [1 -1], b_true, [-1 2 3]), [1 -1 3]), [1 2 -1]);


a_init = randn(P, Q);
b_init = randn(Q, Q, Q);
% x0 = [vec(a_init); vec(b_init)];
% x0 = [vec(a_true); vec(b_true)];

a_est = a_init;
b_est = b_init;

opts = optimset('Algorithm','levenberg-marquardt','Display','on', 'Jacobian', 'on', 'DerivativeCheck', 'off');
for iit = 1:10
  [a_est resnorm] = lsqnonlin(@(x) myfuna(x), vec(a_est), [], [], opts);
  a_est = reshape(a_est, P, Q);

  norm(vec(B-tprod(a_est, [3 -1], tprod(a_est, [2 -1], tprod(a_est, [1 -1], b_est, [-1 2 3]), [1 -1 3]), [1 2 -1])))/norm(vec(B))
  subspace(a_true, a_est)
  
  [b_est resnorm] = lsqnonlin(@(x) myfunb(x), vec(b_est), [], [], opts);
  b_est = reshape(b_est, [Q Q Q]);
  
  norm(vec(B-tprod(a_est, [3 -1], tprod(a_est, [2 -1], tprod(a_est, [1 -1], b_est, [-1 2 3]), [1 -1 3]), [1 2 -1])))/norm(vec(B))
  subspace(a_true, a_est)
  
end

keyboard

function [F J] = myfuna(x)

a = real(reshape(x(1:P*Q), P, Q));
b = b_est;

F = (vec(B-tprod(a, [3 -1], tprod(a, [2 -1], tprod(a, [1 -1], b, [-1 2 3]), [1 -1 3]), [1 2 -1])));

for qq = 1:Q
  aba1(:, :, qq) = a*squeeze(b(qq, :, :))*a';
end
for qq = 1:Q
  aba2(:, :, qq) = a*squeeze(b(:, qq, :))*a';
end
for qq = 1:Q
  aba3(:, :, qq) = a*squeeze(b(:, :, qq))*a';
end

da = zeros(P, P, P, P, Q);
for p1 = 1:P
  for p2 = 1:P
    for p3 = 1:P
      for p4 = 1:P
        for q1 = 1:Q
          if p4 == p1
            da(p1, p2, p3, p4, q1) = da(p1, p2, p3, p4, q1) + aba1(p2, p3, q1);
          end
          if p4 == p2
            da(p1, p2, p3, p4, q1) = da(p1, p2, p3, p4, q1) + aba2(p1, p3, q1);
          end
          if p4 == p3
            da(p1, p2, p3, p4, q1) = da(p1, p2, p3, p4, q1) + aba3(p1, p2, q1);
          end
        end
      end
    end
  end
end

J = -reshape(da, P^3, P*Q);

end



function [F J] = myfunb(x)

a = a_est;
b = reshape(x, [Q Q Q]);

F = (vec(B-tprod(a, [3 -1], tprod(a, [2 -1], tprod(a, [1 -1], b, [-1 2 3]), [1 -1 3]), [1 2 -1])));

J = -kron(kron(a, a), a);

end

end
