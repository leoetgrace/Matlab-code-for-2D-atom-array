%%%This function is used to generate the Green tensor 
%%%Note we need not consider the self-Green function

%%%After the calculation, I find that the size of Green function is related to the distance(r_m)
%%%in the case that if two position is close enough that can be deem as self Green function, it will spoil the final result   



function G=gen_Gt(ka,ri,rj) %ka is the wave vector of the incident field
G=zeros(3); %Green tensor is a 3 by 3 matrix represents different direction.
r=ri-rj; %The relative position of the two atoms ri and rj in 3D space
r_m=sqrt(r(1)^2+r(2)^2+r(3)^2);
f0=exp(1j*ka*r_m)/2/ka/r_m; %The coefficient, note here the coefficient lambda in front of the G is taken into consideration

if abs(r)<10^-3
return
end
for i=1:3
    for k=1:3
        G(i,k)=(-1+(3-1j*3*ka*r_m)/(ka*r_m)^2)*r(i)*r(k)/r_m^2;
        if i==k
            G(i,k)=(1+(1j*ka*r_m-1)/(ka*r_m)^2)+G(i,k);
        end
    end
end
G=f0*G;  
end