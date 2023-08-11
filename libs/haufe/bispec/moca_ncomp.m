function [Fout,wall]=moca_ncomp(F);
% The input F is the set of ns  sources
% F(i,j,k) is the moment of the  k.th source 
% in voxel i in direction j. 
%
% Fout contains the demixed sources. 
% wall is the KxK demixing matrix. 
% 
[ngrid ndum ns]=size(F);
c=zeros(ns,ns);

for i=1:ns;for j=1:ns;
     c(i,j)=(reshape(F(:,:,i),ngrid*3,1))'*reshape(F(:,:,j),ngrid*3,1);
    end;end;



[u s v]=svd(c);
w=u*sqrt(inv(s))*v';



F2=0*F;
for i=1:ns;for j=1:ns;
        F2(:,:,i)=F2(:,:,i)+w(i,j)*F(:,:,j);
end;end

% to check:
% for i=1:ns;for j=1:ns;
%      c2(i,j)=(reshape(F2(:,:,i),ngrid*3,1))'*reshape(F2(:,:,j),ngrid*3,1);
%     end;end;
% c2=c2


F=F2;
nrun=5*ns^2;
wall=w;
for kk=1:nrun
  cost=f2cost(F); 
  
  cmax=0;
  for i=1:ns;
      for j=i+1:ns;
         if cost(i,j)>cmax
             im=i;
             jm=j;
             cmax=cost(i,j);
         end
      end
  end


  
  kkk=0;
  while kkk==0
   im=ceil(rand*ns);
   jm=ceil(rand*ns);
   if im ~= jm
       kkk=1;
   end
  end
  
h1=F(:,:,im);
h2=F(:,:,jm);

M=sum(h1.*h2,2);
N=sum(h2.^2-h1.^2,2)/2;

a=sum(M.^2);
b=sum(M.*N);
c=sum(N.^2);


xmin1=atan(2*b/(a-c))/2;
fmin1=a*cos(xmin1)^2+2*b*cos(xmin1)*sin(xmin1)+c*sin(xmin1)^2;
if xmin1<0
    xmin2=xmin1+pi/2;
else
    xmin2=xmin1-pi/2;
end
fmin2=a*cos(xmin2)^2+2*b*cos(xmin2)*sin(xmin2)+c*sin(xmin2)^2;
if fmin1<fmin2;
    xmin=xmin1;
    fmin=fmin1;
    fmax=fmin2;
else
    xmin=xmin2;
    fmin=fmin2;
    fmax=fmin1;
end
phi=xmin/2;
%disp([im jm phi]) ;
g1=cos(phi)*h1+sin(phi)*h2;
g2=-sin(phi)*h1+cos(phi)*h2;



wall_loc=eye(ns);
wall_loc([im,jm],[im,jm])=[[cos(phi) sin(phi)];[-sin(phi) cos(phi)]];

wall=wall_loc*wall;

F(:,:,im)=g1;
F(:,:,jm)=g2;


 
end

costnew=1000*f2cost(F);
Fout=F;
return;

function cost=f2cost(F) 
[ngrid ndum ns]=size(F);
cost=zeros(ns,ns);
for i=1:ns;for j=1:ns;
        cost(i,j)=sum( (sum( F(:,:,i).*F(:,:,j),2) ).^2); 
end;end;

return

