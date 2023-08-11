function jatjd = calc_jatjd(a,dr)
    [nchan,n] = size(a);
    b = a' * a;
    
    jatjd = zeros(nchan,n,n,n,n);
    
    dloc1 = zeros(n,n,n);
    dloc2 = zeros(n,n,n);
    dloc3 = zeros(n,n,n);
    
    for t = 1:n
        dloc1(t,:,:) = b' * squeeze(dr(t,:,:)) * b;
        dloc2(:,t,:) = b' * squeeze(dr(:,t,:)) * b;
        dloc3(:,:,t) = b' * squeeze(dr(:,:,t)) * b;
    end
    
    for i = 1:n
        for j = 1:n
            for k = 1:n
                jatjd(:,i,:,j,k) = a * dloc1(i,j,k);
            end
        end
    end
    
    for i=1:n
        for j=1:n
            for k=1:n
                jatjd(:,i,j,:,k) = squeeze(jatjd(:,i,j,:,k)) + a * dloc2(j,i,k);
            end
        end
    end
        
    for i=1:n
        for j=1:n
            for k=1:n
                jatjd(:,i,j,k,:) = squeeze(jatjd(:,i,j,k,:)) + a * dloc3(j,k,i);
            end
        end
    end
end