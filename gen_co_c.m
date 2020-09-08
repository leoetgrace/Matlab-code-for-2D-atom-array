%%%This function is used to discover the circular geometry of the atom
%%%array

function co=gen_co_c(N,r)
theta=2*pi/N;
co=zeros(3,N);
for i=1:N
R=r*exp(1j*(i-1)*theta);
co(:,i)=[real(R)+0.5;imag(R)+0.5;0.000001];
end
%co(:,N+1)=[0.5;0.5;0.000001];
end

