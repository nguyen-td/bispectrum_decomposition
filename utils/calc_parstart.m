function [a,d,erstart,model] = calc_parstart(bs,n)
    [nchan,nchan,nchan] = size(bs);
    
    erfit = zeros(3,1); % maybe we should rather loop over n?
    for k = 1:3
        if k == 1
            indspermute = [1,2,3];
        elseif k == 2
            indspermute = [3,1,2];
        else
            indspermute = [2,3,1];
        end
    
        bx = reshape(permute(bs,indspermute),nchan^2,nchan);
        bx = [real(bx); imag(bx)];
        bx2 = bx' * bx;
        [u,s,v] = svd(bx2);
        a = u(:,1:n); % A initialized as the first n left singular vectors
        % a = a_ori;
        % a = a + randn(size(a)) / 10;
        
        jdtjd=calc_jdtjd(a);
        jdtjdr=reshape(jdtjd,n^3,n^3);
        jdtBreal=calc_jdtB(a,real(bs));
        jdtBimag=calc_jdtB(a,imag(bs));
        jdtBrealr=reshape(jdtBreal,n^3,1);
        jdtBimagr=reshape(jdtBimag,n^3,1);
        dreal=inv(jdtjdr)*jdtBrealr;
        dimag=inv(jdtjdr)*jdtBimagr;
        d=dreal+sqrt(-1)*dimag;
        d=reshape(d,n,n,n);
        
        
        scale_a=mean(abs(a(:)));
        scale_b=mean(abs(d(:)));
        lambda=(scale_a/scale_b)^(.25);
        d=d*lambda^3;
        a=a/lambda;
        bs_est=calc_bsmodel(a,d);
        
        
        erfit(k)=norm(bs(:)-bs_est(:))/norm(bs(:));
        model{k}.a=a;
        model{k}.d=d;
        model{k}.er=erfit;
    end
    %erfit=erfit
    [erfitmin,kmin]=min(erfit);
    a=model{kmin}.a;
    d=model{kmin}.d;
    erstart=erfitmin;
end
