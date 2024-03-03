%% Dome with oil 1st Experiment
clc; 
clear all;
close all; 

%% Input Parameters

fontSize = 12;
WaveLength = 0.463e-6;
waveNumber=2*pi/WaveLength;
nDome = 1.46;%1.225

%:----- Dome
NoRing = 8;                             % No. of rings.
NoPhi =  [1 6 12 20 24 32 36 48 48 ]; % No. of LEDs/Ring
thetas = [0 7 14 21 29 37 47 58 73 ]; % Angel for each Ring[deg]

% NoPhi =  [1 6 12 18]; % No. of LEDs/Ring
% thetas = [0 6 12 18]; % Angle for each Ring[deg]


radius = 80e-3;                             % Dome readius mm
ledHeight = -3e-3; %-5.5e-3                         % New sample position
Nled = sum(NoPhi(1:NoRing));                % No. of LEDs

%:----- Camera and Lens
spsize = 2.4e-6/10;                     % CCD pixel size on the object plane including magnification, m
psize = spsize/4;                       % high res pixel size, m,  upscaling ratio chosen
NA = 0.28; % on the object plane 
cutoffFrequency = NA*waveNumber; % rad/m (Not cycle/m!)

%:----- Wave vector
Kx_relative = zeros(length(thetas),max(NoPhi));
Ky_relative = zeros(length(thetas),max(NoPhi));
Kz_relative = zeros(length(thetas),max(NoPhi));

% x y z location of each LED with respect to the dome center

xloc = zeros(length(thetas),max(NoPhi));
yloc = zeros(length(thetas),max(NoPhi));
zloc = zeros(length(thetas),max(NoPhi));

Rotation1 = 0;
Rotation1 = 0; %%47
for j1=1:length(thetas)
    for j2=1:(NoPhi(j1))        
        Phi = thetas(j1);
        Psi = j2*360/NoPhi(j1)
        if j1 > 1
        Psi = j2*360/NoPhi(j1) + Rotation1 ;
        end 
        if j1 > 19
        Psi = j2*360/NoPhi(j1)+ Rotation2 ;
        end
        xloc(j1,j2) =  radius*sind(Phi)*cosd(Psi);
       
        yloc(j1,j2) =  radius*sind(Phi)*sind(Psi);
        zloc(j1,j2) = -radius*cosd(Phi);
    end
end

% Unit Vector and components for each LED to the new sample plane

for j1 = 1:length(thetas)
    for j2=1:(NoPhi(j1))
    
    Kxc = (0 - xloc(j1, j2)); % x componet 
    Kyc = (0 - yloc(j1, j2)); % y component
    Kzc = (ledHeight - zloc(j1, j2));

    V = sqrt(Kxc^2+Kyc^2+Kzc^2); % length of vector  
    Kx_relative(j1, j2) = Kxc/V; % unit vector x component
    Ky_relative(j1, j2) = Kyc/V; % unit vector y component
    Kz_relative(j1, j2) = Kzc/V;
    end
end

kx = nDome*waveNumber*Kx_relative; % full kx of illuminating plane wave
ky = nDome*waveNumber*Ky_relative;

%% Loading images
%:----- Input image size
CenV = 5496/2;
CenH = 3672/2;

%:----- Cropped segment
r = 128*2;
c = 128*2;
I1 = zeros(r,c,Nled);
I2 = zeros(r,c,Nled);


disp("Loading images...");
for j1 = 1:Nled           
    j1
%   Path = strcat('C:\Fourier Ptychography\MyMatlab\Nazabat Algorithms\Dome Algorithm\Dome Dataset 5 Coated Tube Lens\E',int2str(j1),'.tiff');
  Path = strcat('D',int2str(j1),'.tiff');
  SubImageRead = imread(Path);
    I1(:,:,j1) = SubImageRead(CenH - r/2:CenH + r/2-1, CenV - c/2:CenV + c/2-1);
    
% Denoising paramters
alpha = 2;
Beta = 2;
wname = 'db45';
Level = 10;
type = 'h';

  %:----- Noise removing @ 0 Deg
      [C,S] = wavedec2(double(I1(:,:,j1)),Level,wname);
      det1 = Beta*detcoef2('compact',C,S,1);
      sigma1 = median(abs(det1))/0.6745;
      edge1 = wbmpen(C,S,sigma1,alpha);
      I1(:,:,j1) = wthresh(I1(:,:,j1),type,edge1); 
% % % % % -----------
     bk1 = mean2((I1(145:170,164:178,j1)));

    Ibk(j1) = bk1;
    value = 200;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Ibk(j1)>value  % adaptive threshold
        if j1>1
        Ibk(j1) = Ibk(j1-1);
        else
        Ibk(j1) = value;
        end 
    end
    
   I1(:,:,j1)=I1(:,:,j1)-Ibk(j1);
    dd(j1)=min(min(I1(:,:,j1)));
    
end

I1(I1(:,:,:) < 0.00)=0.0; % thresholding and clipping

I1 = sqrt(I1(:,:,:));
disp("Loading done...");

%:----- Correcting Experimental Sequence to Simulation
expSeq =   [
1, ...
7, 2, 3, 4, 5, 6,...
19, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,...
39, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, ...
63, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62,...
95, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, ...
131, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130,...
179, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178,...
227, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226];

for j1 = 1:Nled
    I2(:,:,j1) = I1(:,:,expSeq(j1));
end

%% Square size CCD raw images
m = r*(spsize/psize);
n = c*(spsize/psize);

oxsize=psize*m; % object size, m
oysize=psize*n; % object size, m

oxaxis=-oxsize/2:psize:oxsize/2-psize; % x=0 at matrix location 129
oyaxis=-oysize/2:psize:oysize/2-psize;

[OX,OY]=meshgrid(oxaxis,oyaxis); % x,y axis on object

dkx=2*pi/(psize*(m)); % on the high res. axis plane = 2*kmax/(n) ??? for dome ?
dky=2*pi/(psize*(n)); % 

kmax=pi/spsize; % max rad. freq, for the CCD
dk=2*pi/oxsize;

kx1=-pi/spsize:dkx:pi/spsize-dk; % kx=0 at location 129 for 256x256 matrix.
ky1=-pi/spsize:dky:pi/spsize-dk;
[kxm, kym] = meshgrid(kx1, ky1);
CTF=((kxm.^2+kym.^2) <= cutoffFrequency^2); % the coherent transfer functlon

% figure;imagesc(kx1,ky1,CTF);colorbar;axis image;
% title('CTF');xlabel('kx, rad/m');ylabel('ky, rad/m');

kx2=-pi/psize:dkx:pi/psize-dk; % kx=0 at location 129 for 256x256 matrix.
ky2=-pi/psize:dky:pi/psize-dk;
[kxm2 kym2] = meshgrid(kx2, ky2);
CTFhi=((kxm2.^2+kym2.^2) <= cutoffFrequency^2); % the coherent transfer functlon

%% Synthesized NA and Aliasing 

NA_illum = nDome*(max(max(Kx_relative)))
NA_synth = NA + NA_illum;
spsizemax = WaveLength/2/NA;
psizemax = WaveLength/4/NA_synth;

FPM_res = WaveLength*0.5/NA_synth;

CTF_synth=((kxm2.^2+kym2.^2) <= (NA_synth*waveNumber)^2); % the coherent transfer functlon

%% Image Recovery, basic algorithm
% 
objectRecover=imresize(I2(:,:,1),[m n]); % a possible start up guess
objectRecoverFT=fftshift(fft2(objectRecover));

figure;imagesc(log(abs(objectRecoverFT)));colormap default,colorbar; axis image;
figure;imagesc((abs(objectRecover)));colormap gray,colorbar; axis image;
Pupil=double(CTF);

for loop=1:100
    disp(loop)
    ii=0;
        for tt=1:NoRing
            for nn=1:1:(NoPhi(tt))

               ii=ii+1;
                kxc=round(n/2+1+kx(tt,nn)/dkx); % center of shifted pupil
                kyc=round(m/2+1-ky(tt,nn)/dkx);

                kxl=round(kxc-r/2); % lower and upper limits of shifted pupil
                kxh=round(kxc+r/2)-1;
                kyl=round(kyc-c/2);
                kyh=round(kyc+c/2)-1;

                ObjectRecover_1 = objectRecoverFT(kyl:kyh,kxl:kxh);
                ObjectRecover_2 = ObjectRecover_1.*Pupil;

                ImageLowResolution = ifft2(ifftshift(ObjectRecover_2));
                ImageLowResolution=(m/r)^2*I2(:,:,ii).*exp(1i*angle(ImageLowResolution));

                lowResFT_p=fftshift(fft2(ImageLowResolution));
                
                objectRecoverFT(kyl:kyh,kxl:kxh)=ObjectRecover_1 + ...
                    abs(Pupil).*conj(Pupil.*CTF)./(max(max((abs(Pupil))))).*(lowResFT_p-ObjectRecover_2)./(abs(Pupil).^2+1.0);

                if loop ~= 1.0           
                    Pupil=Pupil + ...
                   0.5*abs(ObjectRecover_1).*conj(ObjectRecover_1)./(max(max((abs(ObjectRecover_1))))).*(lowResFT_p-ObjectRecover_2)./(abs(ObjectRecover_1).^2+1.0);  
                 end

            end
        end

end

objectRecover=ifft2(ifftshift(objectRecoverFT.*CTF_synth)); 

%% Results

figure;
imagesc(log(abs(objectRecoverFT).*CTF_synth));
colormap default,colorbar;
axis image;
title('Fourier Transform of reccovered object', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc((abs(objectRecover))); caxis([0 250])
colormap gray,colorbar;
axis image;
title('Recovered amplitude', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc((angle(objectRecover)));
colormap gray,colorbar; axis image;
title('Recovered phase', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc((angle(Pupil.*CTF)));
colormap default,colorbar; axis image;
title('Recovered Pupil phase', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc((abs(Pupil)));
colormap default,colorbar; axis image;
title('Recovered Pupil ', 'FontSize', fontSize, 'FontName', 'calibri')
pause;
%% tuning of kxc and kyc of each LED
xs=[-30:1:30]; % in pixels, center of shift
ys=[-30:1:30];

ii=0;

        for tt=1:NoRing
            for nn=1:1:(NoPhi(tt))
ii=ii+1
                
    for tz=1:1:length(xs)
    for nz=1:1:length(ys)

                kxc=round(n/2+1+kx(tt,nn)/dkx)+xs(tz); % center of shifted pupil
                kyc=round(m/2+1-ky(tt,nn)/dkx)+ys(nz);

                kxl=round(kxc-r/2); % lower and upper limits of shifted pupil
                kxh=round(kxc+r/2)-1;
                kyl=round(kyc-c/2);
                kyh=round(kyc+c/2)-1;

                ObjectRecover_1 = objectRecoverFT(kyl:kyh,kxl:kxh);
                ObjectRecover_2 = ObjectRecover_1.*Pupil;

                ImageLowResolution = ifft2(ifftshift(ObjectRecover_2));
                I3=abs(ImageLowResolution); % propagated amp          
                I3=I3/max(max(I3)); % normalize
                
                Icap=(m/r)^2*I2(:,:,ii); % captured amp
                Icap=Icap/max(max(Icap)); % normalize
                E2=sum(sum(abs(I3-Icap)))/sum(sum(Icap)); % error function
          
                Err2(tz,nz)=E2;
                Err3(tz,nz)=ssim(I3,Icap);
           
    end
    end

    minMatrix = max(Err3(:));
    [row,col] = find(Err3==minMatrix);

%    minMatrix = min(Err2(:));
%    [row,col] = find(Err2==minMatrix);
    
    kx_LED_cor(tt,nn)=xs(row(1));
    ky_LED_cor(tt,nn)=ys(col(1));

        end
   end

%% recovery with tuned LED positions, work ok, slightly better recovery

%objectRecover=imresize(I2(:,:,1),[m n]); % a possible start up guess to start from zero again
%objectRecoverFT=fftshift(fft2(objectRecover));

%Pupil=double(CTF);

for loop=1:100
    disp(loop)
    ii=0.0;
        for tt=1:NoRing
            for nn=1:1:(NoPhi(tt))

               ii=ii+1;
                kxc=round(n/2+1+kx(tt,nn)/dkx)+kx_LED_cor(tt,nn); % center of shifted pupil
                kyc=round(m/2+1-ky(tt,nn)/dkx)+ky_LED_cor(tt,nn);

                kxl=round(kxc-r/2); % lower and upper limits of shifted pupil
                kxh=round(kxc+r/2)-1;
                kyl=round(kyc-c/2);
                kyh=round(kyc+c/2)-1;

                ObjectRecover_1 = objectRecoverFT(kyl:kyh,kxl:kxh);
                ObjectRecover_2 = ObjectRecover_1.*Pupil;

                ImageLowResolution = ifft2(ifftshift(ObjectRecover_2));
                ImageLowResolution=(m/r)^2*I2(:,:,ii).*exp(1i*angle(ImageLowResolution));

                lowResFT_p=fftshift(fft2(ImageLowResolution));


% Different version of EPRY

                objectRecoverFT(kyl:kyh,kxl:kxh)=ObjectRecover_1 + ...
                    abs(Pupil).*conj(Pupil.*CTF)./(max(max((abs(Pupil))))).*(lowResFT_p-ObjectRecover_2)./(abs(Pupil).^2+1.0);

                if loop ~= 1.0           
                    Pupil=Pupil + ...
                   0.5*abs(ObjectRecover_1).*conj(ObjectRecover_1)./(max(max((abs(ObjectRecover_1))))).*(lowResFT_p-ObjectRecover_2)./(abs(ObjectRecover_1).^2+1.0);  
                end

            end
        end

end

objectRecover=ifft2(ifftshift(objectRecoverFT.*CTF_synth)); 

%% Results

figure;
imagesc(log(abs(objectRecoverFT).*CTF_synth));
colormap default,colorbar;
axis image;
title('Fourier Transform of reccovered object', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc((abs(objectRecover)));
colormap gray,colorbar;
axis image;
title('Recovered amplitude', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc((angle(objectRecover)));
colormap gray,colorbar; axis image;
title('Recovered phase', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc((abs(Pupil)));
colormap default,colorbar; axis image;
title('Recovered pupil phase', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc((angle(Pupil.*CTF)));
colormap default,colorbar; axis image;
title('Recovered pupil', 'FontSize', fontSize, 'FontName', 'calibri')
