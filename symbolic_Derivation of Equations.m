
%%%%%%%%%%%%%%%%%%%%% Symbolic evaluation of Song's derivation

syms theta delta phi Ein Eout E0 I0

Ein=1/sqrt(2)*[1; i]; % LCP light

R=[cos(theta) -sin(theta); sin(theta) cos(theta)] % rotation matrix

P=  [exp(-i*delta/2) 0; 0 exp(+i*delta/2)] % as suggested by Dag

Lphi=[cos(phi)^2 sin(phi)*cos(phi); sin(phi)*cos(phi) sin(phi)^2] % polarization channel on pixels


Js = R*P*inv(R) % Jones matrix of the sample
Js=simplify(Js)


Eout=Lphi*Js*Ein % light on the detector

Iout=Eout(1,1)^2+Eout(2,1)^2 % Intensity on the CMOS
Iout=simplify(Iout)

% phi=0 deg polarization channel, 

L0=subs(Lphi,phi,0)

E0=subs(Eout,phi,0)

I0=abs(E0(1,1))^2+abs(E0(2,1))^2

% Example, delat=0.2 rad, theta= 0.2 rad

I0Song=1/2*(1-sin(.2)*sin(2*0.2)) % as  given by Song equation assuming it is correct

I0mine=1/2*(1+sin(.2)*sin(2*0.2)) % as given by our paper

I0sym=double(subs(I0,[delta theta],[0.2 0.2])) % different from Song ! Hahahahahahahahahahahahah


% Hahhahahahahahahahahahahaahahahhaahahahahahahahahahahahahahahahahahahahahahahahahahah
