function [jtj,jtB]=calc_jtj_and_jtb(a,d,bsdiff)
    [nchan,n] = size(a); 
    jatjareal = calc_jatja(a,real(d));
    jatjarealr = reshape(jatjareal,nchan * n,nchan * n);
    jatjaimag = calc_jatja(a,imag(d)); 
    jatjaimagr = reshape(jatjaimag,nchan * n,nchan * n);
    jatjar = jatjarealr + jatjaimagr;
    
    
    jatjdreal = calc_jatjd(a,real(d));
    jatjdrealr = reshape(jatjdreal,n * nchan,n^3);
    jatjdimag = calc_jatjd(a,imag(d));
    jatjdimagr = reshape(jatjdimag,n * nchan,n^3);
    
    jdtjd = calc_jdtjd(a);
    jdtjdr = reshape(jdtjd,n^3,n^3);
    
    jatBreal = calc_jatB(a,real(d),real(bsdiff));
    jatBrealr = reshape(jatBreal,nchan*n,1);
    
    jdtBreal = calc_jdtB(a,real(bsdiff));
    jdtBrealr = reshape(jdtBreal,n^3,1);
    
    jatBimag = calc_jatB(a,imag(d),imag(bsdiff));
    jatBimagr = reshape(jatBimag,nchan * n,1);
    
    jdtBimag = calc_jdtB(a,imag(bsdiff));
    jdtBimagr = reshape(jdtBimag,n^3,1);
    
    
    jtB = [jatBrealr + jatBimagr; jdtBrealr; jdtBimagr];
    
    jtjda = [jatjdrealr'; jatjdimagr'];
    jtjdd=[[jdtjdr; zeros(size(jdtjdr))], [zeros(size(jdtjdr)); jdtjdr]];
    
    jtj = [[jatjar; jtjda], [jtjda'; jtjdd]];
end