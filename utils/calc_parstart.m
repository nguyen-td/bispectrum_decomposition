function [a,d,erstart,model] = calc_parstart(bs,n)
    [nchan,nchan,nchan,nfreqpairs] = size(bs);
    
    errfit = zeros(3,1);
    for k = 1:3
        if k == 1
            indspermute = [1,2,3,4];
        elseif k == 2
            indspermute = [3,1,2,4];
        else
            indspermute = [2,3,1,4];
        end
    
        bx = reshape(permute(bs,indspermute),nchan^2,nchan,nfreqpairs);
        bx = [real(bx); imag(bx)];
        bx2 = pagemtimes(bx, 'transpose', bx, 'none'); 
        [u,s,v] = pagesvd(bx2);
        a = mean(u(:,1:n,:), 3); % not optimal
        % a = a_ori;
        % a = a + randn(size(a)) / 10;
        
        bs_est = zeros(nchan,nchan,nchan,nfreqpairs);
        d_est = zeros(n,n,n,nfreqpairs);
        for ifreq = 1:nfreqpairs
            jdtjd = calc_jdtjd(a); 
            jdtjdr = reshape(jdtjd,n^3,n^3); 
            jdtBreal = calc_jdtB(a,real(bs(:,:,:,ifreq)));
            jdtBimag = calc_jdtB(a,imag(bs(:,:,:,ifreq)));
            jdtBrealr = reshape(jdtBreal,n^3,1);
            jdtBimagr = reshape(jdtBimag,n^3,1);
            dreal = inv(jdtjdr) * jdtBrealr;
            dimag = inv(jdtjdr)* jdtBimagr;
            d = dreal + 1i * dimag;
            d = reshape(d,n,n,n); 
            
            scale_a = mean(abs(a(:)));
            scale_b = mean(abs(d(:)));
            lambda = (scale_a / scale_b)^(.25);
            d = d * lambda^3;
            a = a / lambda;
            bs_est(:,:,:,ifreq) = calc_bsmodel(a,d);
            d_est(:,:,:,ifreq) = d;
        end
        
        errfit(k) = norm(bs(:) - bs_est(:)) / norm(bs(:));
        model{k}.a = a;
        model{k}.d = d_est;
        model{k}.er = errfit;
    end

    [errfitmin,kmin] = min(errfit);
    a = model{kmin}.a;
    d = model{kmin}.d;
    erstart = errfitmin;
end
