function [d_est,erstart]=calc_parstart_d(bs, a, n, isinit_a)
    [nchan,nchan,nchan,nfreqpairs]=size(bs);
    
    if isinit_a
        [u,s,v]=svd(a);
        a=u(:,1:n);
    end
    
    bs_est = zeros(nchan,nchan,nchan,nfreqpairs);
    d_est = zeros(n,n,n,nfreqpairs);
    for ifreq = 1:nfreqpairs
        jdtjd = calc_jdtjd(a);
        jdtjdr = reshape(jdtjd,n^3,n^3);
        jdtBreal = calc_jdtB(a,real(bs(:,:,:,ifreq)));
        jdtBimag = calc_jdtB(a,imag(bs(:,:,:,ifreq)));
        jdtBrealr = reshape(jdtBreal,n^3,1);
        jdtBimagr = reshape(jdtBimag,n^3,1);
        dreal = inv(jdtjdr)*jdtBrealr;
        dimag = inv(jdtjdr)*jdtBimagr;
        d = dreal + sqrt(-1) * dimag;
        d = reshape(d,n,n,n);
        
        
        scale_a = mean(abs(a(:)));
        scale_b = mean(abs(d(:)));
        lambda = (scale_a/scale_b)^(.25);
        d = d * lambda^3;
        a = a / lambda;
        bs_est(:,:,:,ifreq) = calc_bsmodel(a,d);
        d_est(:,:,:,ifreq) = d;
    end
    
    erstart=norm(bs(:)-bs_est(:))/norm(bs(:));
end
