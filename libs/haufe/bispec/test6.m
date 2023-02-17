load clab_10_20;
load clab_10_10;

clab = clab_10_20;
sa = prepare_sourceanalysis(clab, 'icbm152b_sym');

P = length(clab);
Q_true = 3;
Q_est = 3;

snr = 0.5;
nstart = 1;
i = sqrt(-1);

dips = [sa.grid_coarse(ceil(size(sa.grid_coarse, 1)*rand(Q_true, 1)), :), randn(Q_true, 3)];
a_true = forward_general(dips, sa.fp);
b_true = randn(Q_true, Q_true, Q_true)+i*randn(Q_true, Q_true, Q_true);

B = tprod(a_true, [3 -1], tprod(a_true, [2 -1], tprod(a_true, [1 -1], b_true, [-1 2 3]), [1 -1 3]), [1 2 -1]);
no = randn(size(B)) + i*randn(size(B));
no = norm(vec(B))*no/norm(vec(no));
B = snr*B + (1-snr)*no;

disp('solve with lm using exact Jacobian')
tic
[a_lm b_lm f_lm] = lowrank_bispec(B, Q_est, struct('nstart', nstart, 'solver', 'lm'));
toc
angle_lm = subspace(a_true, a_lm)

disp('solve with lm using exact Jacobian, which is however not explictly stored (memory efficient)')
tic
[a_lm b_lm f_lm] = lowrank_bispec(B, Q_est, struct('nstart', nstart, 'solver', 'lm', 'Jmult', 1));
toc
angle_lm = subspace(a_true, a_lm) 

disp('solve with lm using finite difference approximation')
tic
[a_lmfd b_lmfd f_lmfd] = lowrank_bispec(B, Q_est, struct('nstart', nstart, 'solver', 'lm', 'jacobian', 'off'));
toc
angle_lmd = subspace(a_true, a_lmfd)

disp('solve with l-bfgs')
tic
[a_lbfgs b_lbfgs f_lbfgs] = lowrank_bispec(B, Q_est, struct('nstart', nstart, 'solver', 'lbfgs'));
toc
angle_lbfgs = subspace(a_true, a_lbfgs)
