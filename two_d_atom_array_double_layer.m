%%%This program is used to deploy a numerical solution to the reflection of
%%%a 2D atom array mirror. We deploy a double layer of the atom array

%%%We first start with the trival geometry e.g. square lattice, then we
%%%will go on to discover more about the light in various geometry of the atom array

%%%parameter define
gamma=1; gamma_nr=0; gamma_0=gamma+gamma_nr;%Set the decay rate
delta1=0*gamma; 
delta2=0*gamma; %The detuning between the atom and the incident field
lambda_r=1; %Lambda_a/lambda=1+delta/omega_a 

%%%Geometry set
L=10; N=L*L; %A L by L square lattice
co=gen_co_d(N,0.1,-0.1); %generate the coordinate of the two atom arrays
a_2_l=0.2; %The ratio between a and wavelength 
ka=a_2_l*2*pi; %The wave vector of the incident field

%%%Now we wanna change the angle of the second layer of atoms
psi=0; R=[cos(psi),-sin(psi),0;sin(psi),cos(psi),0;0,0,1]; %Rotation matrix
for i=1:N
    co(:,i+N)=R*co(:,i+N);
end


%%%Green tensor of 3N by 3N matrix
G=zeros(3*(2*N));
for i=1:(2*N)^2
    x=mod((i-1),2*N+1)+1; y=(i-x)/(2*N+1)+1;
    G(3*x-2:3*x,3*y-2:3*y)=gen_Gt(ka,co(:,x),co(:,y));
end

%%%Magic matrix    Also, we wanna introduce atom array with different
%%%detuning
G(:,1:3*N)=3*lambda_r^3/(2*delta1+1j*gamma_0/gamma)*G(:,1:3*N);
G(:,3*N+1:2*3*N)=3*lambda_r^3/(2*delta2+1j*gamma_0/gamma)*G(:,3*N+1:2*3*N); 
I=eye(3*(2*N));
M=I+G; %3*lambda_r^3/(2*delta1+1j*gamma_0/gamma)*G;

%Set E0 to have an angle bewteen the z axis
E0=zeros(3,2*N);  A=1; %Amplitude of the field
w0=0.3*L; %beam waist at its focal point
theta=0;
%theta=pi/3;

x=co(1,:); y=co(2,:); z=co(3,:);
x_c=x*cos(theta)-z*sin(theta); y_c=y; z_c=z*cos(theta)+x*sin(theta);

zr=pi*w0^2;
w=w0*sqrt(1+(z_c*2*pi/ka/zr).^2); %Parameter in the gaussian beam
r=z_c.*(1+(zr./z_c*ka/2/pi).^2);  %%%
phi=atan(z_c*2*pi/ka/zr); 
E0(1,:)=A*w0./w.*exp(1j*ka*z_c).*exp(-1j*phi).*exp(-(x_c.^2+y_c.^2)./w.^2).*exp(1j*ka*(x_c.^2+y_c.^2)/2./r)*cos(theta);
E0(3,:)=-A*w0./w.*exp(1j*ka*z_c).*exp(-1j*phi).*exp(-(x_c.^2+y_c.^2)./w.^2).*exp(1j*ka*(x_c.^2+y_c.^2)/2./r)*sin(theta);
%E0(2,:)=-A*w0./w.*exp(1j*ka*z_c).*exp(-1j*phi).*exp(-(x_c.^2+y_c.^2)/w.^2).*exp(1j*ka*(x_c.^2+y_c.^2)/2./r);
E0=E0(:);

%Take the final step to figure out the field at the atom position
[Q,R]=qr(M);
E=inv(R)*inv(Q)*E0;


%%%plug the E into the input and output equation
%axis range define
Az=3;  Ax=3; %range 0f the x and z plane
A_z=Az/a_2_l;  A_x=Ax/a_2_l;
dd=0.02/a_2_l; %interval of the coordinate 
xx=-A_x:dd:A_x;  zz=-A_z:dd:A_z;  
NN=length(xx)*length(zz); zlen=length(zz); xlen=length(xx);
position=zeros(3,NN);  %space position
for i=1:NN
    x_p=mod(i-1,xlen)+1; z_p=(i-x_p)/xlen+1;
    position(:,i)=[xx(length(xx)-x_p+1);0;zz(z_p)];
end
zz=position(3,:);
zz(zz==0)=0.000001; %set z to be 0.001 in order to avoid the NaN appeared in r
xx=position(1,:); yy=position(2,:);
xx_c=xx*cos(theta)-zz*sin(theta); yy_c=yy;  zz_c=zz*cos(theta)+xx*sin(theta);
%input field
E0_t=zeros(3,NN); 
ww=w0*sqrt(1+(zz_c*2*pi/ka/zr).^2);    %ww is the same
rr=zz_c.*(1+(zr./zz_c*ka/2/pi).^2); %rr is the same
phip=atan(zz_c*2*pi/ka/zr);  %phip is the same

E0_t(1,:)=A*w0./ww.*exp(1j*ka*zz_c).*exp(-1j*phip).*exp(-(xx_c.^2+yy_c.^2)./ww.^2).*exp(1j*ka*(xx_c.^2+yy_c.^2)/2./rr)*cos(theta);
E0_t(3,:)=-A*w0./ww.*exp(1j*ka*zz_c).*exp(-1j*phip).*exp(-(xx_c.^2+yy_c.^2)./ww.^2).*exp(1j*ka*(xx_c.^2+yy_c.^2)/2./rr)*sin(theta);
%E0_t(2,:)=-A*w0./ww.*exp(1j*ka*zz_c).*exp(-1j*phip).*exp(-(xx_c.^2+yy_c.^2)./ww.^2).*exp(1j*ka*(xx_c.^2+yy_c.^2)/2./rr);
%E0_t=E0_t(:); %Input field. note that redundency should be eliminated to improve the calculation
%Multipied by Green tensor which is a NN by 3N matrix

E_f=zeros(3,NN);
parfor i=1:NN
    GG=zeros(3,3*(2*N));
    for j=1:N
        GG(:,3*j-2:3*j)=3*lambda_r^3/(2*delta1+1j*gamma_0/gamma)*gen_Gt(ka,position(:,i),co(:,j));
    end
    for j=N+1:2*N
        GG(:,3*j-2:3*j)=3*lambda_r^3/(2*delta2+1j*gamma_0/gamma)*gen_Gt(ka,position(:,i),co(:,j));
    end
    GG(:,3*(2*N+1)-2:3*(2*N+1))=3*lambda_r^3/(2*delta1+1j*gamma_0/gamma)*gen_Gt(ka,position(:,i),co(:,2*N+1));
    E_f(:,i)=GG*E;
end
E_f=E0_t-E_f;    %%%field at the space

x=position(1,:);
data=find(abs(x-0)<0.000001);
Ex_i=E0_t(1,:); Ex_s=E_f(1,:)-Ex_i; 
%hold on
%plot((-A_z:dd:A_z)*a_2_l,imag(Ex_i(data)))
%plot((-A_z:dd:A_z)*a_2_l,imag(Ex_s(data)))

%E0_a=abs(E0_t(1,:).^2)+abs(E0_t(2,:).^2)+abs(E0_t(3,:).^2);
Ef_a=abs(E_f(1,:)).^2+abs(E_f(2,:)).^2+abs(E_f(3,:)).^2;
comz=abs(E_f(3,:)).^2;
comy=abs(E_f(2,:)).^2;
comx=abs(E_f(1,:)).^2;
%plot((-A_z:dd:A_z)*a_2_l,Ef_a(data))
%T=reshape(I,length(-A_x:dd:A_x),length(-A_z:dd:A_z)); 
T_f=reshape(Ef_a,length(-A_x:dd:A_x),length(-A_z:dd:A_z));
T_f(T_f>4)=0;
imagesc((-A_z:dd:A_z)*a_2_l,(A_x:-dd:-A_x)*a_2_l,T_f)
colorbar

ylabel('x/\lambda')
xlabel('z/\lambda')
%title('L=26, a/\lambda=0.2, \theta=0, b=0.0002')