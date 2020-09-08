%%%This function is used to generate the coordinate od double layers of atoms

function co_r=gen_co_d(N,b1,b2)  %Note that b is the space length between the two array
L=sqrt(N);
co=cell(1,2*L+1);  % we set the lattice constant to be a=1, while ka=2*pi*a/lambda=.....
max=(L-1)/2;
        for i=1:N
            y=mod((i-1),L)+1; x=(i-y)/L+1;
            co{i}=[-max+(x-1),max-(y-1),b1];
        end
        for i=1:N
            y=mod((i-1),L)+1; x=(i-y)/L+1;
            co{i+N}=[-max+(x-1),max-(y-1),b2];
        end
        co_r=zeros(3,2*N);
        for i=1:2*N
            co_r(:,i)=co{i}';
        end
end