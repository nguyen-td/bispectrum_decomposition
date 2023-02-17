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

disp('solve with lm using exact Jacobian')
tic
[a_lm b_lm f_lm] = lowrank_bispec(B, Q, struct('nstart', nstart, 'solver', 'lm'));
toc
angle_lm = subspace(a_true, a_lm)

disp('solve with lm using finite difference approximation')
tic
[a_lmfd b_lmfd f_lmfd] = lowrank_bispec(B, Q, struct('nstart', nstart, 'solver', 'lm', 'jacobian', 'off'));
toc
angle_lmd = subspace(a_true, a_lmfd)

disp('solve with l-bfgs')
tic
[a_lbfgs b_lbfgs f_lbfgs] = lowrank_bispec(B, Q, struct('nstart', nstart, 'solver', 'lbfgs'));
toc
angle_lbfgs = subspace(a_true, a_lbfgs)
