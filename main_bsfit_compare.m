%% load empirical bispectrum bispectrum
load bsori

bs=bsori;
[nchan,nchan,nchan]=size(bs);

n=3; % model order 

%% Stefan 

bx=reshape(bs,nchan^2,nchan);
bx=[real(bx);imag(bx)];
bx2=bx'*bx;
[u,s,v]=svd(bx2);
para=struct('nstart',1, 'solver', 'lm', 'Jmult', 1,'astart',u(:,1:n));
tic
[a_lm b_lm f_lm] = lowrank_bispec_withini(bs,n,para );
toc
bs_est=calc_bsmodel(a_lm,b_lm);
err=norm(bs(:)-bs_est(:))/norm(bs(:))


%% Guido 

tic
[a_lm b_lm err] = bsfit(bs,n);
toc
bs_est=calc_bsmodel(a_lm,b_lm);
err=norm(bs(:)-bs_est(:))/norm(bs(:))
