%%%This fucntion is used to generate the coordinate of the atom array

function co_r=gen_co(N,d)
L=sqrt(N);
co=cell(1,L); a=1; % we set the lattice constant to be a=1, while ka=2*pi*a/lambda=.....
max=(L-1)/2;
for i=1:N
     y=mod((i-1),L)+1; x=(i-y)/L+1;
     if mod(y,2)==0
     co{i}=[-max+(x-1),(max-(y-1)),d];
     else
         co{i}=[(-max+(x-1)),(max-(y-1)),d];
     end
end

 co_r=zeros(3,N); %Reshape the coordinate of the atom
 for i=1:N
    co_r(:,i)=co{i}';
end
 
 