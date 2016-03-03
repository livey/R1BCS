function [t,A,x,flip_loc]=Gen_Data_Flip(M,N,K,K_flip)
x=zeros(N,1);
ind=randperm(N);
ind=ind(1:K);
x(ind)=random('normal',0,1,K,1);
x=x./norm(x);
A=randn(M,N);
for ii=1:N
    A(:,ii)=A(:,ii)/norm(A(:,ii));
end
y=A*x;
t=sgn(y);

ind_flip=randperm(M);
flip_loc=ind_flip(1:K_flip);
t(flip_loc)=~(t(flip_loc));
