n=5; m=5; r=2; J=100; 

rng(2019)

X0=randn(n,r);  Y0=randn(m*J,r); Z=X0*Y0';

XX=randn(n,r*J);  YY=randn(m*J,r); LambdaLambda=randn(n,r*(J-1)); 

rho=100;  eta=1e-4; Ir=repmat(eye(r),J-1,1); 

maxit=1000000;
% Bregman GADMM
Err=zeros(maxit,1); CErr=Err;
X=XX;Y=YY;Lambda=LambdaLambda;
for it=1:maxit
    for j=1:J        
    Xj=X(:,(j-1)*r+1:j*r);
    Yj=Y((j-1)*m+1:j*m,:);
    Zj=Z(:,(j-1)*m+1:j*m);
    tmpy=(norm(Yj(:))^2/2+1)/eta;
    %tmpy=1/eta;
    if j==1  
    X1=(2*(Zj-Xj*Yj')*Yj+(Lambda+rho*(X(:,r+1:end)-repmat(Xj,1,J-1)))*Ir+tmpy*Xj)/tmpy;
    X(:,1:r)=X1;
    else
    Lambdaj=Lambda(n,(j-2)*r+1:(j-1)*r);
    Xj=(2*(Zj-Xj*Yj')*Yj-Lambdaj+rho*(X1-Xj)+tmpy*Xj)/tmpy;
    Lambda(:,(j-2)*r+1:(j-1)*r)=Lambdaj+rho*(Xj-X1);
    X(:,(j-1)*r+1:j*r)=Xj;
    end   
    tmpx=(norm(Xj(:))^2/2+1)/eta;
   % tmpx=1/eta;
    Y((j-1)*m+1:j*m,:)=(2*(Zj-Xj*Yj')'*Xj+tmpx*Yj)/tmpx;
    end   
    
    Errtmp=0; CErrtmp=0; X1=X(:,1:r);
    for j=1:J
        Xj=X(:,(j-1)*r+1:j*r);
        Yj=Y((j-1)*m+1:j*m,:);
        Zj=Z(:,(j-1)*m+1:j*m);
        CErrtmp=CErrtmp+norm(Xj-X1,'fro')^2;
        Errtmp=Errtmp+norm(Zj-Xj*Yj','fro')^2;
    end 
    CErr(it)=CErrtmp;
    Err(it)=Errtmp;    
end
semilogy([CErr Err])
drawnow
    
