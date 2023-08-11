function jatB=calc_jatB(a,dr,Br)
    [nchan,n] = size(a);
    
    jatB = zeros(nchan,n);
    for s = 1:nchan
        for t = 1:n
            jatB(s,t) = ...
            trace(squeeze(dr(t,:,:)) * a' * transpose(squeeze(Br(s,:,:))) * a) ...
        +   trace(squeeze(dr(:,t,:)) * a' * transpose(squeeze(Br(:,s,:))) * a) ... 
        +   trace(squeeze(dr(:,:,t)) * a' * transpose(squeeze(Br(:,:,s))) * a);
        end
    end
end