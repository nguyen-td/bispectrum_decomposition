function [a,d,err,err_all,bsmodel]=bsfit(bs,n,para)
% fits a low rank model to a given cross-spectrum. 
% the model reads:
% bsmodel(p,q,r)=sum_{i,j,k}a(p,i)a(q,j)a(r,k)d(i,j,k) 
% a(p,i) is real valued and d(i,j,k)  complex valued. 
% The indices p,q,r run from 1 to nchan (number of channels)
% and the indices i,j,k run from 1 to n, where n is the model order. 
% Of course, n is lower than nchan. 
% The parameters a can be considered as a mixing matrix, and d is the cross-bispectrum 
% of the n sources. Note, that any nxn transformation of the mixing matrix  can be absorbed into d, 
% and hence only the space spanned by the mixing matrix, and not the matrix  itself,  is unique.  
%  
%
% Input:
% bs: nchan x nchan x nchan in general complex tensor cross-spectral values for nchan channels 
% n:  model order
% para: (optional) 
% para.a (starting values for the mixing matrix a, e.g.
%         para.a=randn(nchan,n);
%         If no starting value is given, it is estimated from bs
%
% Output
% a: nchan x n mixing matrix
% d: source cross-bispectrum is size nxnxn
% err: relative model error, i.e. norm(bs-bsmodel)/norm(bs) with norm
%        meaning Frobenius norm
% errall: errors for all iteration steps (just to check the progress)
% bsmodel: the model  cross-bispectrum 


% defaults: 
alpha=.01; % starting value for regularation of LM procedure
kmax=50; %maximum number of iterations 
kmin=8; % minimum number of iterations 
a=[];  %starting value for mixing matrix. If empty, it is estimated from bs
errtol=10^(-6); % programs stops if error decreases by less than this value 
                % within the last two iterations. Only calculated after 7
                % iterations. (For bad starting values 
               

[nchan,nchan,nchan]=size(bs);
if nargin>2;
    if isfield(para,'a')
        a=para.a;
    end
end

if length(a)==0
  [a,d,erstart]=calc_parstart(bs,n);
else
    [d,erstart]=calc_parstart_d(bs,a,n);
end


err=erstart;
err_all=zeros(kmax+1,1);
err_all(1)=err;
 
kont=1;
k=0;
while kont==1;
    k=k+1;
    bs_est=calc_bsmodel(a,d);
    bsdiff=bs-bs_est;
    [jtj,jtB]=calc_jtj_and_jtb(a,d,bsdiff);
    npar=length(jtj);
    jtj_regu=jtj+alpha*trace(jtj)*eye(npar)/npar;
    par_new=inv(jtj_regu)*jtB;
    a_new=a+reshape(par_new(1:nchan*n),nchan,n);
    d_new_real=real(d)+reshape(par_new(nchan*n+1:nchan*n+n^3),n,n,n);
    d_new_imag=imag(d)+reshape(par_new(nchan*n+1+n^3:end),n,n,n);
    d_new=d_new_real+sqrt(-1)*d_new_imag;
    bs_est_new=calc_bsmodel(a_new,d_new);
     err_new=norm(bs(:)-bs_est_new(:))/norm(bs(:));
     if err_new<err
         alpha=alpha/10;
         alpha=max([10^(-8),alpha]);
         a=a_new;
         d=d_new;
         err=err_new;
     else
         alpha=alpha*10;
     end
    % log(alpha)/log(10)
     err_all(k+1)=err;
     if k==kmax;kont=0;end
     if k>kmin-1;
         diff_err=err_all(k-1)-err_all(k+1);
         if diff_err<errtol;kont=0;end
     end
     bs_est=calc_bsmodel(a,d);

end
 
err_all=err_all(1:k+1);

bsmodel=bs_est;
return;

function [jtj,jtB]=calc_jtj_and_jtb(a,d,bsdiff)

[nchan,n]=size(a);
jatjareal=calc_jatja(a,real(d));
jatjarealr=reshape(jatjareal,nchan*n,nchan*n);
jatjaimag=calc_jatja(a,imag(d));
jatjaimagr=reshape(jatjaimag,nchan*n,nchan*n);
jatjar=jatjarealr+jatjaimagr;


jatjdreal=calc_jatjd(a,real(d));
jatjdrealr=reshape(jatjdreal,n*nchan,n^3);
jatjdimag=calc_jatjd(a,imag(d));
jatjdimagr=reshape(jatjdimag,n*nchan,n^3);

jdtjd=calc_jdtjd(a);
jdtjdr=reshape(jdtjd,n^3,n^3);

jatBreal=calc_jatB(a,real(d),real(bsdiff));
jatBrealr=reshape(jatBreal,nchan*n,1);

jdtBreal=calc_jdtB(a,real(bsdiff));
jdtBrealr=reshape(jdtBreal,n^3,1);

jatBimag=calc_jatB(a,imag(d),imag(bsdiff));
jatBimagr=reshape(jatBimag,nchan*n,1);

jdtBimag=calc_jdtB(a,imag(bsdiff));
jdtBimagr=reshape(jdtBimag,n^3,1);


jtB=[jatBrealr+jatBimagr;jdtBrealr;jdtBimagr];

jtjda=[jatjdrealr';jatjdimagr'];
jtjdd=[[jdtjdr;zeros(size(jdtjdr))],[zeros(size(jdtjdr));jdtjdr]];

jtj=[[jatjar;jtjda],[jtjda';jtjdd]];

return;

function jatB=calc_jatB(a,dr,Br)

[nchan,n]=size(a);

jatB=zeros(nchan,n);
for s=1:nchan;for t=1:n;
        jatB(s,t)=...
        trace(squeeze(dr(t,:,:))*a'*transpose(squeeze(Br(s,:,:)))*a)...
    +   trace(squeeze(dr(:,t,:))*a'*transpose(squeeze(Br(:,s,:)))*a)... 
    +   trace(squeeze(dr(:,:,t))*a'*transpose(squeeze(Br(:,:,s)))*a);
    end;end


return;

function jatja=calc_jatja(a,dr)
[nchan,n]=size(a);
b=a'*a;


E1=zeros(n,n,n,n);
E2=zeros(n,n,n,n);
E3=zeros(n,n,n,n);

for i=1:n;for k=1:n;
    E1(i,:,k,:)=transpose(squeeze(dr(:,i,:)))*b*squeeze(dr(:,k,:));
    E2(i,:,k,:)=transpose(squeeze(dr(i,:,:)))*b*squeeze(dr(k,:,:));
    E3(i,:,k,:)=squeeze(dr(i,:,:))*b*transpose(squeeze(dr(k,:,:)));
 end;end
  

   

   jatja=zeros(nchan,n,nchan,n);
   
  
   for t=1:n;for w=1:n
           jatja(:,t,:,w)=squeeze(jatja(:,t,:,w))...
               +(a*squeeze(E3(t,:,:,w))*transpose(a))' ...
               +(a*squeeze(E2(t,:,:,w))*transpose(a))'...
               + (a*squeeze(E3(:,t,w,:))*transpose(a))'...
               + (a*squeeze(E1(t,:,:,w))*transpose(a))'...
               + (a*squeeze(E2(:,t,w,:))*transpose(a))'...
               + (a*squeeze(E1(:,t,w,:))*transpose(a))';
       end;end
   
   E1cont=zeros(n,n);
   E2cont=zeros(n,n);
   E3cont=zeros(n,n);
   for t=1:n;for w=1:n;
           E1cont(t,w)=trace(squeeze(E1(t,:,w,:))*b');
           E2cont(t,w)=trace(squeeze(E2(:,t,:,w))*b');
           E3cont(t,w)=trace(squeeze(E3(t,:,w,:))*b');
    end;end

   for s=1:nchan;
        jatja(s,:,s,:)=squeeze(jatja(s,:,s,:))...
            +E1cont+E2cont+E3cont;
   end
   
   return;
   
   function jatjd=calc_jatjd(a,dr)

[nchan,n]=size(a);
b=a'*a;

jatjd=zeros(nchan,n,n,n,n);

dloc1=zeros(n,n,n);
dloc2=zeros(n,n,n);
dloc3=zeros(n,n,n);

for t=1:n;
    dloc1(t,:,:)=b'*squeeze(dr(t,:,:))*b;
    dloc2(:,t,:)=b'*squeeze(dr(:,t,:))*b;
    dloc3(:,:,t)=b'*squeeze(dr(:,:,t))*b;
end

for i=1:n;for j=1:n;for k=1:n; 
      jatjd(:,i,:,j,k)=a*dloc1(i,j,k);
end;end;end

for i=1:n;for j=1:n;for k=1:n; 
      jatjd(:,i,j,:,k)=squeeze(jatjd(:,i,j,:,k))+a*dloc2(j,i,k);
end;end;end
    
for i=1:n;for j=1:n;for k=1:n; 
      jatjd(:,i,j,k,:)=squeeze(jatjd(:,i,j,k,:))+a*dloc3(j,k,i);
end;end;end

return;

function JdtB=calc_jdtB(a,B)
[nchan,n]=size(a);
X=zeros(n,n,nchan);
JdtB=zeros(n,n,n);
for i=1:nchan;
    X(:,:,i)=a'*B(:,:,i)*a;
end
if n>1
for i=1:n;
    JdtB(i,:,:)=squeeze(X(i,:,:))*a;
end
else
    JdtB(1,:,:)=transpose(squeeze(X(1,:,:)))*a;
end

    
%JdtB=reshape(JdtB,n^3,1);
return;

function jdtjd=calc_jdtjd(a)

[nchan,n]=size(a);
b=a'*a;



bb=reshape(kron(b,b),n,n,n,n);
t2=zeros(n,n,n,n,n,n);
for p1=1:n;
       for p2=1:n;
           t2(p1,:,:,p2,:,:)=b(p1,p2)*bb;
        end
end

%jdtjd=reshape(t2,n^3,n^3);
jdtjd=t2;
return;

function [a,d,erstart,model]=calc_parstart(bs,n)

[nchan,nchan,nchan]=size(bs);

erfit=zeros(3,1);
for k=1:3;
    if k==1;
        indspermute=[1,2,3];
    elseif k==2;
        indspermute=[3,1,2];
    else
        indspermute=[2,3,1];
    end
bx=reshape(permute(bs,indspermute),nchan^2,nchan);
bx=[real(bx);imag(bx)];
bx2=bx'*bx;
[u,s,v]=svd(bx2);
a=u(:,1:n);
%a=a_ori;
%a=a+randn(size(a))/10;



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

return;



function [d,erstart]=calc_parstart_d(bs,a,n)

[nchan,nchan,nchan]=size(bs);


[u,s,v]=svd(a);
a=u(:,1:n);


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


erstart=norm(bs(:)-bs_est(:))/norm(bs(:));


return;

function bs=calc_bsmodel(a,d)

[nchan,n]=size(a);
Bx=zeros(n,nchan,nchan);
for i=1:n;
    Bx(i,:,:)=a*squeeze(d(i,:,:))*a';
end
bs=reshape(a*reshape(Bx,n,nchan^2),nchan,nchan,nchan);


return;

