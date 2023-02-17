function [a b resnorm] = lowrank_bispec_withini(B, Q, pars)

% Stefan Haufe, 2013

i = sqrt(-1);
P = size(B, 1);

nstart = 1;
solver = 'lm';
Jmult = 0;
algorithm = 'levenberg-marquardt';
jacobmult = [];
jacobmult_b = [];
sqrtcov = eye(P);
kin=0;
if nargin > 2
  if isfield(pars, 'nstart')
    nstart = pars.nstart;
  end
  if isfield(pars, 'solver')
    solver = pars.solver;
  end
  if isequal(solver, 'lm') && isfield(pars, 'Jmult')
    Jmult = pars.Jmult;
  end
  if isfield(pars, 'sqrtcov')
    sqrtcov = pars.sqrtcov;
  end
   if isfield(pars, 'astart')
    astart = pars.astart;
    kin=1;
  end
end
  
if Jmult
  jacobmult = @myfun_Jmult;
  jacobmult_b = @myfun_Jmult_b;
  algorithm = 'trust-region-reflective';
end

sqrtcov = sqrtcov / norm(sqrtcov);

nob = norm(vec(B));
B = B / nob;

a = zeros(P, Q);
b = zeros(Q^3);
aba1 = zeros(P, P, Q);
aba2 = zeros(P, P, Q);
aba3 = zeros(P, P, Q);

% @optimplotresnorm
opts = optimset('Algorithm', algorithm, 'PlotFcns', [], 'Display','off', 'Jacobian', 'on', 'FinDiffType', 'central', 'DerivativeCheck', 'off', 'JacobMult', jacobmult_b);
%   opts = optimset('Display','off', 'GradObj', 'on', 'FinDiffType', 'central', 'DerivativeCheck', 'off');

% tic

for istart = 1:nstart*10
     a_init(:, istart) = vec(orth(sqrtcov*randn(P, Q)));
 
    if kin==1
            if istart==1
                 a_init(:, istart)=vec(astart);
              
            end
    end
     a_tmp = vec(a_init(:,istart));
  [b_init(:, istart), resnorm(istart)] = lsqnonlin(@(x) myfun_lm_b_real(x), zeros(2*Q^3, 1), [], [], opts);
%   [b_init(:, istart), resnorm(istart)] = fminunc(@(x) myfun_lbfgs_b(x), zeros(2*Q^3, 1), opts);
end

resnorm=resnorm

% toc

[so, in] = sort(resnorm);
a_init = a_init(:, in(1:nstart));
b_init = b_init(:, in(1:nstart));
resnorm = resnorm(in(1:nstart));  

if isequal(solver, 'lm')
  % levenberg-marquardt
  % trust-region-reflective
  % 'PlotFcns', @optimplotresnorm, 
  opts = optimset('Algorithm', algorithm, 'PlotFcns', [], 'Display','off', 'Jacobian', 'on', 'FinDiffType', 'central', 'DerivativeCheck', 'off', 'JacobMult', jacobmult); 

%   tic
  for istart = 1:nstart
%     disp(istart)
% tic
    [x(:, istart) resnorm(istart)] = lsqnonlin(@(x) myfun_lm_real(x), [a_init(:, istart); b_init(:, istart)], [], [], opts);
% toc
  end
%   toc
else
  % l-BFGS
  opts = optimset('Display','off', 'GradObj', 'on', 'FinDiffType', 'central', 'DerivativeCheck', 'off');

  % tic
  for istart = 1:nstart
%     disp(istart)
    [x(:, istart) resnorm(istart)] = fminunc(@(x) myfun_lbfgs(x), [a_init(:, istart); b_init(:, istart)], opts);
  end
  % toc
  
end

[f in] = min(resnorm);
x = x(:, in);
a = reshape(x(1:P*Q), P, Q);
b = reshape(x((P*Q+1):end), [Q Q Q 2]);
b = b(:, :, :, 1) + i*b(:, :, :, 2);

b = b * nob;

% computes unsquared errors and Jacobian for levenberg marquard
function [f J] = myfun_lm(x)

  a = real(reshape(x(1:P*Q), P, Q));
  b = reshape(x((P*Q+1):end), [Q Q Q]);

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
  
  J = [reshape(J, P^3, P*Q) -kron(kron(a, a), a)];
  
end

% computes unsquared errors and Jacobian for levenberg marquard
function [f J] = myfun_lm_real(x)

%   disp(1)
  a = reshape(x(1:P*Q), P, Q);
  b = reshape(x((P*Q+1):end), [Q Q Q 2]);
  b = b(:, :, :, 1) + i*b(:, :, :, 2);

  f = vec(B-tprod(a, [3 -1], tprod(a, [2 -1], tprod(a, [1 -1], b, [-1 2 3]), [1 -1 3]), [1 2 -1]));
  f = [real(f); imag(f)];    
  
  for qq = 1:Q
    aba1(:, :, qq) = a*squeeze(b(qq, :, :))*a';
    aba2(:, :, qq) = a*squeeze(b(:, qq, :))*a';
    aba3(:, :, qq) = a*squeeze(b(:, :, qq))*a';
  end
    
  J = sparse(2*P^3, P*Q + 2*Q^3);

  if ~Jmult
    Ja = zeros(P, P, P, P, Q);
    for pp = 1:P
      Ja(pp, :, :, pp, :) = -aba1;
      Ja(:, pp, :, pp, :) = squeeze(Ja(:, pp, :, pp, :)) - aba2;
      Ja(:, :, pp, pp, :) = squeeze(Ja(:, :, pp, pp, :)) - aba3;
    end

    Ja = reshape(Ja, P^3, P*Q);
    Jb = -kron(kron(a, a), a);
    
    J = sparse([real(Ja) Jb zeros(size(Jb));  imag(Ja) zeros(size(Jb)) Jb]);
  end
   
end



function JY = myfun_Jmult(Jinfo, Y, flag)

  siY = size(Y);

  if flag >= 0
  %     JY = J*Y;
    JY = zeros(2*P^3, siY(2));
    for pp = 1:P
      for qq = 1:Q
        J = zeros(P, P, P);
        J(pp, :, :) = -aba1(:, :, qq);
        J(:, pp, :) = squeeze(J(:, pp, :)) - aba2(:, :, qq);
        J(:, :, pp) = squeeze(J(:, :, pp)) - aba3(:, :, qq);
        JY = JY + [vec(real(J)); vec(imag(J))]*Y(sub2ind([P Q], pp, qq), :);
      end
    end
    for q1 = 1:Q
      for q2 = 1:Q
        for q3 = 1:Q
          kro = -kron(kron(a(:, q3), a(:, q2)), a(:, q1));
          JY = JY + [kro; sparse(P^3, 1)]*Y(P*Q+sub2ind([Q Q Q], q1, q2, q3), :);
          JY = JY + [sparse(P^3, 1); kro]*Y(P*Q+Q^3+sub2ind([Q Q Q], q1, q2, q3), :);
        end
      end
    end
  end
  
  if flag == 0
%     JY = J'*J*Y;
    Y = JY;
  end
  
  siY = size(Y);
  
  if flag <= 0
%     JY = J'*Y;
    JYa = zeros(siY(2), P, Q);
    for pp = 1:P
      for qq = 1:Q
        J = zeros(P, P, P);
        J(pp, :, :) = -aba1(:, :, qq);
        J(:, pp, :) = squeeze(J(:, pp, :)) - aba2(:, :, qq);
        J(:, :, pp) = squeeze(J(:, :, pp)) - aba3(:, :, qq);
        JYa(:, pp, qq) = Y'*[vec(real(J)); vec(imag(J))];
      end
    end
    JYbre = zeros(siY(2), Q, Q, Q);
    JYbim = zeros(siY(2), Q, Q, Q);
    for q1 = 1:Q
      for q2 = 1:Q
        for q3 = 1:Q
          kro = -kron(kron(a(:, q3), a(:, q2)), a(:, q1));
          JYbre(:, q1, q2, q3) = Y'*[kro; sparse(P^3, 1)];
          JYbim(:, q1, q2, q3) = Y'*[sparse(P^3, 1); kro];
        end
      end
    end
    JY = [reshape(JYa, siY(2), []) reshape(JYbre, siY(2), []) reshape(JYbim, siY(2), [])]';
  end
  
end

% computes summed squared errors and gradient for lBFGS algorithm
% memory efficient implementation
function [f grad] = myfun_lbfgs(x)

  a = reshape(x(1:P*Q), P, Q);
  b = reshape(x((P*Q+1):end), [Q Q Q 2]);
  b = b(:, :, :, 1) + i*b(:, :, :, 2);

  f = vec(B-tprod(a, [3 -1], tprod(a, [2 -1], tprod(a, [1 -1], b, [-1 2 3]), [1 -1 3]), [1 2 -1]));

  J = zeros(P, P, P);
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



% computes unsquared errors and Jacobian for levenberg marquard for fixed a
function [f J] = myfun_lm_b(x)

  a = reshape(a_tmp, [P Q]);
  b = reshape(x, [Q Q Q]);

  f = (vec(B-tprod(a, [3 -1], tprod(a, [2 -1], tprod(a, [1 -1], b, [-1 2 3]), [1 -1 3]), [1 2 -1])));

  J = -kron(kron(a, a), a);
  
end

% computes unsquared errors and Jacobian for levenberg marquard for fixed a
function [f J] = myfun_lm_b_real(x)

  a = reshape(a_tmp, [P Q]);
  b = reshape(x, [Q Q Q 2]);
  b = b(:, :, :, 1) + i*b(:, :, :, 2);

  f = vec(B-tprod(a, [3 -1], tprod(a, [2 -1], tprod(a, [1 -1], b, [-1 2 3]), [1 -1 3]), [1 2 -1]));
  f = [real(f); imag(f)];   
  
  J = sparse(2*P^3, 2*Q^3);
  if ~Jmult
    Jb = -kron(kron(a, a), a);
    J = [Jb zeros(size(Jb)); zeros(size(Jb)) Jb];
  end
end

function JY = myfun_Jmult_b(Jinfo, Y, flag)
  siY = size(Y);
  
  if flag >= 0
  %     JY = J*Y;
    JY = zeros(2*P^3, siY(2));
    for q1 = 1:Q
      for q2 = 1:Q
        for q3 = 1:Q
          kro = -kron(kron(a(:, q3), a(:, q2)), a(:, q1));
          JY = JY + [kro; sparse(P^3, 1)]*Y(sub2ind([Q Q Q], q1, q2, q3), :);
          JY = JY + [sparse(P^3, 1); kro]*Y(Q^3+sub2ind([Q Q Q], q1, q2, q3), :);
        end
      end
    end
  end
  
  if flag == 0
%     JY = J'*J*Y;
    Y = JY;
  end
  
  siY = size(Y);
  
  if flag <= 0
%     JY = J'*Y;
    JYbre = zeros(siY(2), Q, Q, Q);
    JYbim = zeros(siY(2), Q, Q, Q);
    for q1 = 1:Q
      for q2 = 1:Q
        for q3 = 1:Q
          kro = -kron(kron(a(:, q3), a(:, q2)), a(:, q1));
          JYbre(:, q1, q2, q3) = Y'*[kro; sparse(P^3, 1)];
          JYbim(:, q1, q2, q3) = Y'*[sparse(P^3, 1); kro];
        end
      end
    end
    JY = [reshape(JYbre, siY(2), []) reshape(JYbim, siY(2), [])]';
  end
  
end


% computes summed squared errors and gradient for lBFGS algorithm for fixed
% a --  memory efficient implementation
function [f grad] = myfun_lbfgs_b(x)

  a = reshape(a_tmp, P, Q);
  b = reshape(x, [Q Q Q 2]);
  b = b(:, :, :, 1) + i*b(:, :, :, 2);

  f = vec(B-tprod(a, [3 -1], tprod(a, [2 -1], tprod(a, [1 -1], b, [-1 2 3]), [1 -1 3]), [1 2 -1]));
  
  gradb = zeros(Q, Q, Q);
  for q1 = 1:Q
    for q2 = 1:Q
      for q3 = 1:Q
        gradb(q1, q2, q3) = 2*f'*kron(kron(a(:, q3), a(:, q2)), a(:, q1));
      end
    end
  end
  gradb = vec(gradb);
  
  grad = [-real(gradb); imag(gradb)];

  f = sum(real(f).^2 + imag(f).^2);    
end


function v = vec(x)
  v = reshape(x, [], 1);  
end

end
