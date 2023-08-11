function JdtB=calc_jdtB(a,B)
    [nchan,n] = size(a);
    X = zeros(n,n,nchan);
    JdtB = zeros(n,n,n);

    for i = 1:nchan
        X(:,:,i) = a' * B(:,:,i) * a;
    end

    if n > 1
        for i = 1:n
            JdtB(i,:,:) = squeeze(X(i,:,:)) * a;
        end
    else
        JdtB(1,:,:) = transpose(squeeze(X(1,:,:))) * a;
    end
end