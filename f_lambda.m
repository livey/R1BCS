function y=f_lambda(x)
N=length(x);
y=zeros(N,1);
for ii=1:N
   if x(ii)
       y(ii)=1/4/x(ii)*tanh(x(ii)/2);
   else                 %% if x(ii) == 0 , use the infinite approximation 
       y(ii)=1/8;    
   end
end
