function bs=calc_bsmodel(a,d)

[nchan,n]=size(a);
Bx=zeros(n,nchan,nchan);
for i=1:n;
    Bx(i,:,:)=a*squeeze(d(i,:,:))*a';
end
bs=reshape(a*reshape(Bx,n,nchan^2),nchan,nchan,nchan);


return;



