function jdtjd = calc_jdtjd(a)
    [nchan,n] = size(a);
    b = a' * a;
    
    bb = reshape(kron(b,b),n,n,n,n);
    t2 = zeros(n,n,n,n,n,n);
    for p1 = 1:n
        for p2=1:n
            t2(p1,:,:,p2,:,:) = b(p1,p2) * bb; 
        end
    end
    
    jdtjd = t2;
end

