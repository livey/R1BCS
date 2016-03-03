function [x_est,mun] = R1BCS(t,A,max_iter,tao)
% Init 
[M,N]=size(A);
par = 1e-10;
a=par;
b=par;
c=par;
d=par;
% Ealpha=ones(N,1)*1e-3;
Ealpha=ones(N,1);
% Ebeta =ones(M,1)*1e-3;
Ebeta =ones(M,1);
epsilong = ones(M,1);   %% need adjust 
mun = zeros(M,1);
sigman= eye(M);
iter =0;
while (iter < max_iter )
   iter = iter + 1 ;
   
   % Update sigmax mux 
   De=f_lambda(epsilong);
   sigmax = inv(diag(Ealpha)+2*A'*diag(De)*A);
   mux = sigmax*A'*(.5*(2*t-1)-2*diag(De)*mun-2*tao*De);
   
   % Update mun sigman 
   sigman = inv(diag(Ebeta)+2*diag(De));
   mun = sigman*(.5*(2*t-1)-2*diag(De)*A*mux-2*tao*De);
   
%    mun=zeros(M,1);
%    sigman=zeros(M,M);
%    
   % Update Ealpha 
   at = a+.5;
   bt=b+.5*(mux.*mux+diag(sigmax));
   Ealpha = at./bt;
   
   % Update Ebeta 
   ct= c+ .5;
   dt=d+.5*(mun.*mun+diag(sigman));
   Ebeta = ct./dt;
   
   % Max Step , Update epsilong 
   B = mux*mux'+ sigmax ;
   epsilong = diag(A*B*A')+mun.*mun+...
       diag(sigman) + tao*tao*ones(M,1)+...
       2*diag(mun)*A*mux + tao*A*mux + 2*tao*mun;
   epsilong=sqrt(epsilong);
end
x_est=mux/norm(mux);
mun=mun/norm(mun);



