function jatja = calc_jatja(a,dr)
    [nchan,n] = size(a);
    b = a' * a;
    
    E1 = zeros(n,n,n,n);
    E2 = zeros(n,n,n,n);
    E3 = zeros(n,n,n,n);
    
    for i = 1:n
        for k=1:n
            E1(i,:,k,:) = transpose(squeeze(dr(:,i,:))) * b * squeeze(dr(:,k,:));
            E2(i,:,k,:) = transpose(squeeze(dr(i,:,:))) * b * squeeze(dr(k,:,:));
            E3(i,:,k,:) = squeeze(dr(i,:,:)) * b * transpose(squeeze(dr(k,:,:)));
        end
    end
    jatja = zeros(nchan,n,nchan,n);
    
    for t = 1:n
        for w = 1:n
            jatja(:,t,:,w) = squeeze(jatja(:,t,:,w))...
                + (a * squeeze(E3(t,:,:,w)) * transpose(a))' ...
                + (a * squeeze(E2(t,:,:,w)) * transpose(a))'...
                + (a * squeeze(E3(:,t,w,:)) * transpose(a))'...
                + (a * squeeze(E1(t,:,:,w)) * transpose(a))'...
                + (a * squeeze(E2(:,t,w,:)) * transpose(a))'...
                + (a * squeeze(E1(:,t,w,:)) * transpose(a))';
        end
    end
    
    E1cont = zeros(n,n);
    E2cont = zeros(n,n);
    E3cont = zeros(n,n);
    for t = 1:n
        for w = 1:n
            E1cont(t,w) = trace(squeeze(E1(t,:,w,:)) * b');
            E2cont(t,w) = trace(squeeze(E2(:,t,:,w)) * b');
            E3cont(t,w) = trace(squeeze(E3(t,:,w,:)) * b');
        end
    end
    
    for s = 1:nchan
        jatja(s,:,s,:) = squeeze(jatja(s,:,s,:))...
            + E1cont + E2cont + E3cont;
    end
end