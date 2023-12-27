% fits a low rank model to a given cross-spectrum. 
% the model reads:
% bsmodel(p,q,r)=sum_{i,j,k}a(p,i)a(q,j)a(r,k)d(i,j,k) 
% a(p,i) is real valued and d(i,j,k)  complex valued. 
% The indices p,q,r run from 1 to nchan (number of channels)
% and the indices i,j,k run from 1 to n, where n is the model order. 
% Of course, n is lower than nchan. 
% The parameters a can be considered as a mixing matrix, and d is the cross-bispectrum 
% of the n sources. Note, that any nxn transformation of the mixing matrix  can be absorbed into d, 
% and hence only the space spanned by the mixing matrix, and not the matrix  itself,  is unique.  
%  
%
% Input:
% bs: nchan x nchan x nchan in general complex tensor cross-spectral values for nchan channels 
% n:  model order
% para: (optional) 
% para.a (starting values for the mixing matrix a, e.g.
%         para.a=randn(nchan,n);
%         If no starting value is given, it is estimated from bs
% para.isfit_a (boolean, 1 if the mixing matrix a should be fitted, 0 if 
%         the mixing matrix should be held fixed and only d should be fitted 
%         which turns the problem into a convex one)
%
% Output
% a: nchan x n mixing matrix
% d: source cross-bispectrum is size nxnxn
% err: relative model error, i.e. norm(bs-bsmodel)/norm(bs) with norm
%        meaning Frobenius norm
% errall: errors for all iteration steps (just to check the progress)
% bsmodel: the model  cross-bispectrum 

function [a,d,err,err_all,bsmodel] = bsfit(bs,n,para)
    % defaults
    alpha = .01; % starting value for regularation of LM procedure
    kmax = 50;   % maximum number of iterations 
    kmin = 8;    % minimum number of iterations 
    a = [ ];     % starting value for mixing matrix. If empty, it is estimated from bs
    errtol = 10^(-6); % programs stops if error decreases by less than this value 
                      % within the last two iterations. Only calculated after 
                      % 7 iterations. (For bad starting values) 
                   
    [nchan,nchan,nchan] = size(bs);
    isfit_a = true; % fit a by default unless explicitely specified otherwise
    if nargin > 2
        if isfield(para,'a')
            a = para.a;
        end
        if isfield(para, 'isfit_a')
            isfit_a = para.isfit_a;
        end
    end
    
    % initialization
    if isempty(a)
        [a,d,erstart] = calc_parstart(bs,n);
    else
        [d,erstart] = calc_parstart_d(bs, a, n, isfit_a);
    end
    
    err = erstart;
    err_all = zeros(kmax+1,1);
    err_all(1) = err;
     
    kont = 1;
    k = 0;
    while kont == 1
        k = k+1;
        bs_est = calc_bsmodel(a,d);
        bsdiff = bs - bs_est;
    
        [jtj,jtB] = calc_jtj_and_jtb(a,d,bsdiff);
        npar = length(jtj);
        jtj_regu = jtj + alpha * trace(jtj) * eye(npar) / npar;
        par_new = inv(jtj_regu) * jtB;
        
        if isfit_a
            a_new = a + reshape(par_new(1:nchan * n),nchan,n);
        else
            a_new = a;
        end
        d_new_real = real(d) + reshape(par_new(nchan * n + 1:nchan * n + n^3),n,n,n);
        d_new_imag = imag(d) + reshape(par_new(nchan * n + 1 + n^3:end),n,n,n);
        d_new = d_new_real + 1i * d_new_imag;
        bs_est_new = calc_bsmodel(a_new,d_new);
        err_new = norm(bs(:) - bs_est_new(:)) / norm(bs(:));
        if err_new < err
            alpha = alpha/10;
            alpha = max([10^(-8),alpha]);
            a = a_new;
            d = d_new;
            err = err_new;
        else
            alpha = alpha*10;
        end
        % log(alpha)/log(10)
        err_all(k+1) = err;
        if k == kmax
            kont = 0;
        end
        if k > kmin - 1
            diff_err = err_all(k-1) - err_all(k+1);
            if diff_err < errtol
                kont = 0;
            end
        end
        bs_est = calc_bsmodel(a,d);
    end
    err_all = err_all(1:k+1);
    bsmodel = bs_est;
end
   


