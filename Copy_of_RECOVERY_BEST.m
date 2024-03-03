%% Dome 1st Experiment
% Fourier Ptychography algorithm written by Mahdieh
% Department of Microsystem and Maritime Communication 
% University Southeast of Norway - Vestfold Campus
% Date 10/22/2021
% Version 1.0

clc; 
clear;
close all; 

%% Input Parameters

fontSize = 8;
waveLength = 0.525e-6;
WaveNumber = 2*pi/waveLength;
nDome = 1;%1.45;

%:----- Dome
NoRing =  9;%7                             % No. of rings.
%NoPhi =  [1 6 12 18 24 30 36 42 48 54]; % No. of LEDs/Ring
%thetas = [0 6 12 18 25 32 40 50 62 75]; % Angel for each Ring[deg]

NoPhi =  [1 6 12 18 24 30 36 42 48]% 54 30 36 30 24]% 30]; % No. of LEDs/Ring
thetas = [0 6 12 18 25 32 40 50 62]%75 32 4032 25]% 32]; % Angle for each Ring[deg]

radius = 120e-3;                             % Dome readius mm
LedHeight = 16.6e-3%16.6e-3;                           % New sample position
Nled = sum(NoPhi(1:NoRing));                 % No. of LEDs

%:----- Camera and Lens
spsize = 3.45e-6/10;                     % CCD pixel size on the object plane including magnification, m
psize = spsize/4;                       % high res pixel size, m,  upscaling ratio chosen
NA = 0.28; % on the object plane 
CutoffFrequency = NA*WaveNumber; % rad/m (Not cycle/m!)

%%
%:----- Wave vector
Kx_relative = zeros(length(thetas),max(NoPhi));
Ky_relative = zeros(length(thetas),max(NoPhi));
Kz_relative = zeros(length(thetas),max(NoPhi));

% x y z location of each LED with respect to the dome center

Xloc = zeros(length(thetas),max(NoPhi));
Yloc = zeros(length(thetas),max(NoPhi));
Zloc = zeros(length(thetas),max(NoPhi));

Rotation = -5.9;
for j1=1:length(thetas)
    for j2=1:(NoPhi(j1))        
        Phi = thetas(j1);
        if j1 > 1
        Psi = j2*360/NoPhi(j1) + Rotation ;
        else
        Psi = j2*360/NoPhi(j1) ;
        end
        
        Xloc(j1,j2) =  radius*sind(Phi)*cosd(Psi);
        Yloc(j1,j2) =  radius*sind(Phi)*sind(Psi);
        Zloc(j1,j2) = -radius*cosd(Phi);
    end
end

% Unit Vector and components for each LED to the new sample plane

for j1 = 1:length(thetas)
    for j2=1:(NoPhi(j1))
    
    Kxc = (0 - Xloc(j1, j2)); % x componet 
    Kyc = (0 - Yloc(j1, j2)); % y component
    Kzc = (LedHeight - Zloc(j1, j2));

    V = sqrt(Kxc^2+Kyc^2+Kzc^2); % length of vector  
    Kx_relative(j1, j2) = Kxc/V; % unit vector x component
    Ky_relative(j1, j2) = Kyc/V; % unit vector y component
    Kz_relative(j1, j2) = Kzc/V;
    end
end

kx = nDome*WaveNumber*Kx_relative; % full kx,ky of illuminating plane wave
ky = nDome*WaveNumber*Ky_relative;

%% Loading images
%:----- Input image size
% CenV = 5496/2;
% CenH = 3672/2;
CenV = 2448/2;
CenH = 2048/2;

%:----- Cropped segment
r = 256/2;
c = 256/2;
IHH = zeros(r,c,Nled);
IHV = zeros(r,c,Nled);
IVH = zeros(r,c,Nled);
IVV = zeros(r,c,Nled);

imSeqLowRes1 = zeros(r,c,Nled);
imSeqLowRes2 = zeros(r,c,Nled);
imSeqLowRes3 = zeros(r,c,Nled);
imSeqLowRes4 = zeros(r,c,Nled);


imSeqLowRes10 = zeros(r,c,Nled);
imSeqLowRes20 = zeros(r,c,Nled);
imSeqLowRes30 = zeros(r,c,Nled);
imSeqLowRes40 = zeros(r,c,Nled);


Iccd = zeros(2048,2448);
Isub = zeros(2048/2,2448/2);
I0 = zeros(2048,2448);
I45 = zeros(2048,2448);
I90 = zeros(2048,2448);
I135 = zeros(2048,2448);


disp("Loading images...");

for j1 = 1:Nled           
    j1
 
  %:------ Image cropping

  Path = strcat(int2str(j1),'.tiff');
  Iccd = double(imread(Path));

  %:----- Reading the images
  
% MATLAB code for interp2()
% % Specifying a 2-D grid
% [X,Y] = meshgrid(-4:4);
%  
% % Sampling the peaks() function
% V = peaks(X,Y);
% [Xq,Yq] = meshgrid(-4:0.25:4);
%  
% % Calling the interp2() function
% title('Linear Interpolation Using Finer Grid');
% Vq = interp2(X,Y,V,Xq,Yq);

Isub(1:1:end,1:1:end)=Iccd(1:2:end,1:2:end);
I0(1:end-1,1:end-1)=interp2(Isub(:,:),'cubic'); 

Isub(1:1:end,1:1:end)=Iccd(1:2:end,2:2:end);
I45(1:end-1,1:end-1)=interp2(Isub,'cubic');

Isub(1:1:end,1:1:end)=Iccd(2:2:end,2:2:end);
I90(1:end-1,1:end-1)=interp2(Isub,'cubic');

Isub(1:1:end,1:1:end)=Iccd(2:2:end,1:2:end);
I135(1:end-1,1:end-1)=interp2(Isub,'cubic');

  IHH(:,:,j1) = (I0(CenH - r/2 :CenH + r/2-1 , CenV - c/2 :CenV + c/2-1 ));
  IHV(:,:,j1) = (I45(CenH - r/2 :CenH + r/2-1 , CenV - c/2 :CenV + c/2-1 ));
  IVH(:,:,j1) = (I90(CenH - r/2 :CenH + r/2-1 , CenV - c/2 :CenV + c/2-1 ));
  IVV(:,:,j1) = (I135(CenH - r/2 :CenH + r/2-1 , CenV - c/2 :CenV + c/2-1 ));
 
% 
% Denoising paramters
alpha = 2;
Beta = 8;
wname = 'db45';
Level = 10;
type = 'h';

  %:----- Noise removing @ 0 Deg
      [C,S] = wavedec2(double(IHH(:,:,j1)),Level,wname);
      det1 = Beta*detcoef2('compact',C,S,1);
      sigma1 = median(abs(det1))/0.6745;
      edge1 = wbmpen(C,S,sigma1,alpha);
      IHH(:,:,j1) = wthresh(IHH(:,:,j1),type,edge1); 
    %:----- Noise removing @ 45 Deg
      [C,S] = wavedec2(double(IHV(:,:,j1)),Level,wname);
      det1 = Beta*detcoef2('compact',C,S,1);
      sigma1 = median(abs(det1))/0.6745;
      edge1 = wbmpen(C,S,sigma1,alpha);
      IHV(:,:,j1) = wthresh(IHV(:,:,j1),type,edge1); 
    %:----- Noise removing @ 90 Deg
      [C,S] = wavedec2(double(IVH(:,:,j1)),Level,wname);
      det1 = Beta*detcoef2('compact',C,S,1);
      sigma1 = median(abs(det1))/0.6745;
      edge1 = wbmpen(C,S,sigma1,alpha);
      IVH(:,:,j1) = wthresh(IVH(:,:,j1),type,edge1); 
    %:----- Noise removing @ 135 Deg
      [C,S] = wavedec2(double(IVV(:,:,j1)),Level,wname);
      det1 = Beta*detcoef2('compact',C,S,1);
      sigma1 = median(abs(det1))/0.6745;
      edge1 = wbmpen(C,S,sigma1,alpha);
      IVV(:,:,j1) = wthresh(IVV(:,:,j1),type,edge1); 
      
      % Creating adaptive threshold
      
%      bk1HH = mean2((IHH(110:128,60:80,j1)));
%       IbkHH(j1) = bk1HH;
%       %:-----
%       value =120;
%       if IbkHH(j1)>value  % adaptive threshold
%         if j1>1
%         IbkHH(j1) = IbkHH(j1-1);
%         else
%         IbkHH(j1) = value;
%         end
%       end
%       %-----
%       bk1HV = mean2((IHV(110:120,60:80,j1)));
%       IbkHV(j1) = bk1HV;
%       %:-----
%       if IbkHV(j1)>value  % adaptive threshold
%         if j1>1
%         IbkHV(j1) = IbkHV(j1-1);
%         else
%         IbkHV(j1) = value;
%         end
%       end
%       %-----
%       bk1VH = mean2((IVH(110:120,60:80,j1)));
%       IbkVH(j1) = bk1VH;
%       %:-----
%       if IbkVH(j1)>value  % adaptive threshold
%         if j1>1
%         IbkVH(j1) = IbkVH(j1-1);
%         else
%         IbkVH(j1) = value;
%         end
%       end
%       %-----
%       bk1VV = mean2((IVV(110:120,60:80,j1)));
%       IbkVV(j1) = bk1VV;
%       %:-----
%       if IbkVV(j1)>value  % adaptive threshold
%         if j1>1
%         IbkVV(j1) = IbkVV(j1-1);
%         else
%         IbkVV(j1) = value;
%         end
%       end
%       %-----
%       
%       IHH(:,:,j1)=IHH(:,:,j1)-IbkHH(j1);
%       IHV(:,:,j1)=IHV(:,:,j1)-IbkHV(j1);
%       IVH(:,:,j1)=IVH(:,:,j1)-IbkVH(j1);
%       IVV(:,:,j1)=IVV(:,:,j1)-IbkVV(j1);
             
end

%:---- Removing negative values
IHH(IHH(:,:,:) < 0.00) = 0.0;
IHV(IHV(:,:,:) < 0.00) = 0.0;
IVH(IVH(:,:,:) < 0.00) = 0.0;
IVV(IVV(:,:,:) < 0.00) = 0.0;

%:----- Intensity to Amplitude
imSeqLowRes10 = sqrt(IHH(:,:,:));
imSeqLowRes20 = sqrt(IHV(:,:,:));
imSeqLowRes30 = sqrt(IVH(:,:,:));
imSeqLowRes40 = sqrt(IVV(:,:,:));

disp("Loading done...");


%:----- Correcting Experimental Sequence to Simulation
expSeq =   [1,...
            7, 2, 3, 4, 5, 6,...
            19, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,...
            37, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36,...
            61,	38,	39,	40,	41,	42,	43,	44,	45,	46,	47,	48,	49,	50,	51,	52,	53,	54,	55,	56,	57,	58,	59,	60,...
            91,	62,	63,	64,	65,	66,	67,	68,	69,	70,	71,	72,	73,	74,	75,	76,	77,	78,	79,	80,	81,	82,	83,	84,	85,	86,	87,	88,	89,	90,...
            127, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126,...
            169, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168,...
            217, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216,...
            271, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270];

for j1 = 1:Nled
    imSeqLowRes1(:,:,j1) = imSeqLowRes10(:,:,expSeq(j1));
    imSeqLowRes2(:,:,j1) = imSeqLowRes20(:,:,expSeq(j1));
    imSeqLowRes3(:,:,j1) = imSeqLowRes30(:,:,expSeq(j1));
    imSeqLowRes4(:,:,j1) = imSeqLowRes40(:,:,expSeq(j1));
    
end

clear imSeqLowRes10 imSeqLowRes20 imSeqLowRes30 imSeqLowRes40 IHH IHV IVH IVV

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

kx1 = -pi/spsize:dkx:pi/spsize-dk; % kx=0 at location 129 for 256x256 matrix.
ky1 = -pi/spsize:dky:pi/spsize-dk;
[kxm, kym] = meshgrid(kx1, ky1);
CTF = ((kxm.^2+kym.^2) <= CutoffFrequency^2); % the coherent transfer functlon
Pupil = CTF;

% figure;imagesc(kx1,ky1,CTF);colorbar;axis image;
% title('CTF');xlabel('kx, rad/m');ylabel('ky, rad/m');

kx2 = -pi/psize:dkx:pi/psize-dk; % kx=0 at location 129 for 256x256 matrix.
ky2 = -pi/psize:dky:pi/psize-dk;
[kxm2 kym2] = meshgrid(kx2, ky2);
%CTFhi = ((kxm2.^2+kym2.^2) <= CutoffFrequency^2); % the coherent transfer functlon

%% Synthesized NA and Aliasing 

NA_illum = nDome*(max(max(Kx_relative)));
NA_synth = NA + NA_illum;
spsizemax = waveLength/2/NA;
psizemax = waveLength/4/NA_synth;

FPM_res = waveLength*0.5/NA_synth;

CTF_synth=((kxm2.^2+kym2.^2) <= (NA_synth*WaveNumber)^2); % the coherent transfer functlon

%% Image Recovery, Channel 1
%:-----  
objectRecoverHH = imresize(imSeqLowRes1(:,:,1),[m n]); 
objectRecoverFTHH = fftshift(fft2(objectRecoverHH));
figure;imagesc((abs(objectRecoverHH)));
colormap gray;colorbar;
axis image;
title('Low-Res Image,0 deg', 'FontSize', fontSize, 'FontName', 'calibri')


loop = 250;
pupil = double(CTF);
for ji=1:loop
    disp(ji)
    j2 = 0;
    for j3=1:NoRing
            for j4=1:1:(NoPhi(j3))
 
            j2 = j2 + 1;
            %:----- Center of shifted pupil
            kxc=round(n/2+1-kx(j3,j4)/dkx); 
            kyc=round(m/2+1-ky(j3,j4)/dkx);
            
            %:----- Lower and upper limits of shifted pupil
            kxl = round(kxc-r/2); 
            kxh = round(kxc+r/2) - 1;
            kyl = round(kyc-c/2);
            kyh = round(kyc+c/2) - 1;

            ObjectRecover_1HH = objectRecoverFTHH(kyl:kyh,kxl:kxh);
            ObjectRecover_2HH = ObjectRecover_1HH.*pupil.*CTF;
            ImageLowResolutionHH = ifft2(ifftshift(ObjectRecover_2HH));
            ImageLowResolutionHH=(m/r)^2*imSeqLowRes1(:,:,j2).*exp(1i*angle(ImageLowResolutionHH));
            lowResFT_pHH = fftshift(fft2(ImageLowResolutionHH));

% Zheng            
            objectRecoverFTHH(kyl:kyh,kxl:kxh)=objectRecoverFTHH(kyl:kyh,kxl:kxh) + conj(pupil.*CTF)./(max(max((abs(pupil)).^2))).*(lowResFT_pHH - ObjectRecover_2HH);
 
 %           objectRecoverFTHH(kyl:kyh,kxl:kxh)=ObjectRecover_1HH + ...
 %               abs(pupil).*conj(pupil.*CTF)./(max(max((abs(pupil))))).*(lowResFT_pHH - ObjectRecover_2HH)./(abs(pupil).^2+1.0);

            tmp= objectRecoverFTHH(kyl:kyh,kxl:kxh);
            
            if ji ~= 1.0  
           %    pupil = pupil + conj(ObjectRecover_1HH)./(max(max((abs(ObjectRecover_1HH)).^2))).*(lowResFT_pHH - ObjectRecover_2HH);       % Better ?       
          %     pupil = pupil + conj(ObjectRecover_2HH)./(max(max((abs(ObjectRecover_2HH)).^2))).*(lowResFT_pHH - ObjectRecover_2HH);        % better?     
                   pupil = pupil + conj(tmp)./(max(max((abs(tmp)).^2))).*(lowResFT_pHH - ObjectRecover_2HH);        % better     
               

            end 
        end
    end
end

objectRecoverHH = ifft2(ifftshift(objectRecoverFTHH));

figure;
imagesc((abs(objectRecoverHH)));
colormap gray;colorbar;
axis image;
title('Recovered amplitude,0 deg', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc((angle(objectRecoverHH)));
colormap gray;colorbar; axis image;
title('Recovered phase,0 deg', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc(log(abs(objectRecoverFTHH).*CTF_synth));axis image;
title('Fourier Transform of reccovered object 0 deg', 'FontSize', fontSize, 'FontName', 'calibri')
colormap default;colorbar; 

figure;
imagesc((abs(pupil.*CTF)));colorbar; axis image;

figure;
imagesc((angle(pupil.*CTF)));colorbar; axis image;
 pause
%% tuning of kxc and kyc of each LED, does not help

xs=[-5:1:5]; % in pixels, center of pupil shift for each LED
ys=[-5:1:5];

ii=0;

        for tt=1:NoRing
            for nn=1:1:(NoPhi(tt))
ii=ii+1
                
    for tz=1:1:length(xs)
    for nz=1:1:length(ys)

                kxc=round(n/2+1-kx(tt,nn)/dkx)+xs(tz); % center of shifted pupil
                kyc=round(m/2+1-ky(tt,nn)/dkx)+ys(nz);

                kxl=round(kxc-r/2); % lower and upper limits of shifted pupil
                kxh=round(kxc+r/2)-1;
                kyl=round(kyc-c/2);
                kyh=round(kyc+c/2)-1;

                ObjectRecover_1 = objectRecoverFTHH(kyl:kyh,kxl:kxh);
                ObjectRecover_2 = ObjectRecover_1.*pupil.*CTF;

                ImageLowResolution = ifft2(ifftshift(ObjectRecover_2));
                I3=abs(ImageLowResolution); % propagated amp          
                I3=I3/max(max(I3)); % normalize
                
                Icap=(m/r)^2*imSeqLowRes1(:,:,ii); % captured amp
                Icap=Icap/max(max(Icap)); % normalize
                E2=sum(sum(abs(I3-Icap)))/sum(sum(Icap)); % error function
          
                Err2(tz,nz)=E2;
                Err3(tz,nz)=ssim(I3,Icap);
           
            end
    end

         minMatrix = max(Err3(:));
        [row,col] = find(Err3==minMatrix);

        minMatrix = min(Err2(:));
        [row,col] = find(Err2==minMatrix);
    
        kx_LED_cor(tt,nn)=xs(row);
        ky_LED_cor(tt,nn)=ys(col);

            end
       end

%% Recovery with tuned LED positions, does not help

loop = 250;

for ji=1:loop
    disp(ji)
    j2 = 0;
    for j3=1:NoRing
            for j4=1:1:(NoPhi(j3))
 
            j2 = j2 + 1;
            %:----- Center of shifted pupil
            kxc=round(n/2+1-kx(j3,j4)/dkx)+kx_LED_cor(j3,j4); 
            kyc=round(m/2+1-ky(j3,j4)/dkx)+ky_LED_cor(j3,j4);
            
            %:----- Lower and upper limits of shifted pupil
            kxl = round(kxc-r/2); 
            kxh = round(kxc+r/2) - 1;
            kyl = round(kyc-c/2);
            kyh = round(kyc+c/2) - 1;

            ObjectRecover_1HH = objectRecoverFTHH(kyl:kyh,kxl:kxh);
            ObjectRecover_2HH = ObjectRecover_1HH.*pupil.*CTF;
            ImageLowResolutionHH = ifft2(ifftshift(ObjectRecover_2HH));
            ImageLowResolutionHH=(m/r)^2*imSeqLowRes1(:,:,j2).*exp(1i*angle(ImageLowResolutionHH));
            lowResFT_pHH = fftshift(fft2(ImageLowResolutionHH));
            
            objectRecoverFTHH(kyl:kyh,kxl:kxh)=objectRecoverFTHH(kyl:kyh,kxl:kxh) + conj(pupil.*CTF)./(max(max((abs(pupil)).^2))).*(lowResFT_pHH - ObjectRecover_2HH);
 
        end
    end
end
objectRecoverHH = ifft2(ifftshift(objectRecoverFTHH)); 

figure;
imagesc((abs(objectRecoverHH)));
colormap gray;colorbar;
axis image;
title('Recovered amplitude,0 deg', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc((angle(objectRecoverHH)));
colormap gray;colorbar; axis image;
title('Recovered phase,0 deg', 'FontSize', fontSize, 'FontName', 'calibri')


%% Image Recovery, Channel 2
%:-----  
objectRecoverHV = imresize(imSeqLowRes2(:,:,1),[m n]); 
objectRecoverFTHV = fftshift(fft2(objectRecoverHV));
figure;imagesc((abs(objectRecoverHV)));
colormap gray;colorbar;
axis image;
title('Low-Res Image,45 deg', 'FontSize', fontSize, 'FontName', 'calibri')


pupil = double(CTF);
for ji=1:loop
    disp(ji)
    j2 = 0;
    for j3=1:NoRing
            for j4=1:1:(NoPhi(j3))
 
            j2 = j2 + 1;
            %:----- Center of shifted pupil
            kxc=round(n/2+1-kx(j3,j4)/dkx); 
            kyc=round(m/2+1-ky(j3,j4)/dkx);
            
            %:----- Lower and upper limits of shifted pupil
            kxl = round(kxc-r/2); 
            kxh = round(kxc+r/2) - 1;
            kyl = round(kyc-c/2);
            kyh = round(kyc+c/2) - 1;

            ObjectRecover_1HV = objectRecoverFTHV(kyl:kyh,kxl:kxh);
            ObjectRecover_2HV = ObjectRecover_1HV.*pupil.*CTF;
            ImageLowResolutionHV = ifft2(ifftshift(ObjectRecover_2HV));
            ImageLowResolutionHV=(m/r)^2*imSeqLowRes2(:,:,j2).*exp(1i*angle(ImageLowResolutionHV));
            lowResFT_pHV = fftshift(fft2(ImageLowResolutionHV));
            
            objectRecoverFTHV(kyl:kyh,kxl:kxh)=objectRecoverFTHV(kyl:kyh,kxl:kxh) + conj(pupil.*CTF)./(max(max((abs(pupil)).^2))).*(lowResFT_pHV - ObjectRecover_2HV);

            tmp= objectRecoverFTHV(kyl:kyh,kxl:kxh);

            if ji ~= 1.0  
%                pupil = pupil + conj(ObjectRecover_1HV)./(max(max((abs(ObjectRecover_1HV)).^2))).*(lowResFT_pHV - ObjectRecover_2HV);              
                pupil = pupil + conj(tmp)./(max(max((abs(tmp)).^2))).*(lowResFT_pHV - ObjectRecover_2HV);        % better     

            end 
        end
    end
end

objectRecoverHV = ifft2(ifftshift(objectRecoverFTHV));

figure;
imagesc((abs(objectRecoverHV)));
colormap gray;colorbar;
axis image;
title('Recovered amplitude,0 deg', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc((angle(objectRecoverHV)));
colormap gray;colorbar; axis image;
title('Recovered phase,0 deg', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc(log(abs(objectRecoverFTHV).*CTF_synth));axis image;
title('Fourier Transform of reccovered object 0 deg', 'FontSize', fontSize, 'FontName', 'calibri')
colormap default;colorbar; 

figure;
imagesc((abs(pupil.*CTF)));colorbar; axis image;

figure;
imagesc((angle(pupil.*CTF)));colorbar; axis image;

%% Image Recovery, Channel 3
%:-----  
objectRecoverVH = imresize(imSeqLowRes3(:,:,1),[m n]); 
objectRecoverFTVH = fftshift(fft2(objectRecoverVH));
figure;imagesc((abs(objectRecoverVH)));
colormap gray;colorbar;
axis image;
title('Low-Res Image,90 deg', 'FontSize', fontSize, 'FontName', 'calibri')


for ji=1:loop
    disp(ji)
    j2 = 0;
    for j3=1:NoRing
            for j4=1:1:(NoPhi(j3))
 
            j2 = j2 + 1;
            %:----- Center of shifted pupil
            kxc=round(n/2+1-kx(j3,j4)/dkx); 
            kyc=round(m/2+1-ky(j3,j4)/dkx);
            
            %:----- Lower and upper limits of shifted pupil
            kxl = round(kxc-r/2); 
            kxh = round(kxc+r/2) - 1;
            kyl = round(kyc-c/2);
            kyh = round(kyc+c/2) - 1;

            ObjectRecover_1VH = objectRecoverFTVH(kyl:kyh,kxl:kxh);
            ObjectRecover_2VH = ObjectRecover_1VH.*pupil.*CTF;
            ImageLowResolutionVH = ifft2(ifftshift(ObjectRecover_2VH));
            ImageLowResolutionVH=(m/r)^2*imSeqLowRes3(:,:,j2).*exp(1i*angle(ImageLowResolutionVH));
            lowResFT_pVH = fftshift(fft2(ImageLowResolutionVH));
            
            objectRecoverFTVH(kyl:kyh,kxl:kxh)=objectRecoverFTVH(kyl:kyh,kxl:kxh) + conj(pupil.*CTF)./(max(max((abs(pupil)).^2))).*(lowResFT_pVH - ObjectRecover_2VH);
            
            tmp = objectRecoverFTVH(kyl:kyh,kxl:kxh);
            if ji ~= 1.0  
%                pupil = pupil + conj(ObjectRecover_1VH)./(max(max((abs(ObjectRecover_1VH)).^2))).*(lowResFT_pVH - ObjectRecover_2VH);              
                pupil = pupil + conj(tmp)./(max(max((abs(tmp)).^2))).*(lowResFT_pVH - ObjectRecover_2VH);        % better     
 %               
            end 
        end
    end
end
objectRecoverVH = ifft2(ifftshift(objectRecoverFTVH)); 

figure;
imagesc((abs(objectRecoverVH)));
colormap gray;colorbar;
axis image;
title('Recovered amplitude,0 deg', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc((angle(objectRecoverVH)));
colormap gray;colorbar; axis image;
title('Recovered phase,0 deg', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc(log(abs(objectRecoverFTVH).*CTF_synth));axis image;
title('Fourier Transform of reccovered object 0 deg', 'FontSize', fontSize, 'FontName', 'calibri')
colormap default;colorbar; 

figure;
imagesc((abs(pupil.*CTF)));colorbar; axis image;

figure;
imagesc((angle(pupil.*CTF)));colorbar; axis image;


%% Image Recovery, Channel 4
%:----- 
objectRecoverVV = imresize(imSeqLowRes4(:,:,1),[m n]); 
objectRecoverFTVV = fftshift(fft2(objectRecoverVV));
figure;
imagesc((abs(objectRecoverVV)));
colormap gray;colorbar;
axis image;
title('Low-Res Image,135 deg', 'FontSize', fontSize, 'FontName', 'calibri')


for ji=1:loop
    disp(ji)
    j2 = 0;
    for j3=1:NoRing
            for j4=1:1:(NoPhi(j3))
 
            j2 = j2 + 1;
            %:----- Center of shifted pupil
            kxc=round(n/2+1-kx(j3,j4)/dkx); 
            kyc=round(m/2+1-ky(j3,j4)/dkx);
            
            %:----- Lower and upper limits of shifted pupil
            kxl = round(kxc-r/2); 
            kxh = round(kxc+r/2) - 1;
            kyl = round(kyc-c/2);
            kyh = round(kyc+c/2) - 1;

            ObjectRecover_1VV = objectRecoverFTVV(kyl:kyh,kxl:kxh);
            ObjectRecover_2VV = ObjectRecover_1VV.*pupil.*CTF;
            ImageLowResolutionVV = ifft2(ifftshift(ObjectRecover_2VV));
            ImageLowResolutionVV=(m/r)^2*imSeqLowRes4(:,:,j2).*exp(1i*angle(ImageLowResolutionVV));
            lowResFT_pVV = fftshift(fft2(ImageLowResolutionVV));
            
            objectRecoverFTVV(kyl:kyh,kxl:kxh)=objectRecoverFTVV(kyl:kyh,kxl:kxh) + conj(pupil.*CTF)./(max(max((abs(pupil)).^2))).*(lowResFT_pVV - ObjectRecover_2VV);
            
            tmp= objectRecoverFTVV(kyl:kyh,kxl:kxh);
            
            if ji ~= 1.0  
%                pupil = pupil + conj(ObjectRecover_1VV)./(max(max((abs(ObjectRecover_1VV)).^2))).*(lowResFT_pVV - ObjectRecover_2VV);              
                pupil = pupil + conj(tmp)./(max(max((abs(tmp)).^2))).*(lowResFT_pVV - ObjectRecover_2VV);        % better     
 
            end 
        end
    end
end
objectRecoverVV = ifft2(ifftshift(objectRecoverFTVV)); 

figure;
imagesc((abs(objectRecoverVV)));
colormap gray,colorbar;
axis image;
title('Recovered amplitude,135 deg', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc((angle(objectRecoverVV)));
colormap gray,colorbar; axis image;
title('Recovered phase,135 deg', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc(log(abs(objectRecoverFTVV).*CTF_synth));axis image;
title('Fourier Transform of reccovered object 135 deg', 'FontSize', fontSize, 'FontName', 'calibri')
colormap default,colorbar; 

figure;
imagesc((abs(pupil.*CTF)));colorbar; axis image;

figure;
imagesc((angle(pupil.*CTF)));colorbar; axis image;


%% Plot the results

figure;
imagesc((abs(objectRecoverHH)));
colormap gray,colorbar;
axis image;
title('Recovered amplitude,45 deg', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc((angle(objectRecoverHH)));
colormap gray,colorbar; axis image;
title('Recovered phase,45 deg', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc((abs(objectRecoverHV)));
colormap gray,colorbar;
axis image;
title('Recovered amplitude,45 deg', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc((angle(objectRecoverHV)));
colormap gray,colorbar; axis image;
title('Recovered phase,45 deg', 'FontSize', fontSize, 'FontName', 'calibri')


figure;
imagesc((abs(objectRecoverVH)));
colormap gray,colorbar;
axis image;
title('Recovered amplitudex,90 deg', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc((angle(objectRecoverVH)));
colormap gray,colorbar; axis image;
title('Recovered phase,90 deg', 'FontSize', fontSize, 'FontName', 'calibri')


figure;
imagesc((abs(objectRecoverVV)));
colormap gray,colorbar;
axis image;
title('Recovered amplitude, 135 deg', 'FontSize', fontSize, 'FontName', 'calibri')

figure;
imagesc((angle(objectRecoverVV)));
colormap gray,colorbar; axis image;
title('Recovered phase, 135 deg', 'FontSize', fontSize, 'FontName', 'calibri')
colormap default,colorbar;

figure;
imagesc(log(abs(objectRecoverFTHH).*CTF_synth)); axis image
title('Fourier Transform of reccovered object 0 deg', 'FontSize', fontSize, 'FontName', 'calibri')
colormap default,colorbar;
figure;
imagesc(log(abs(objectRecoverFTHV).*CTF_synth)); axis image
title('Fourier Transform of reccovered object 45 deg', 'FontSize', fontSize, 'FontName', 'calibri')
colormap default,colorbar;
figure;
imagesc(log(abs(objectRecoverFTVH).*CTF_synth)); axis image
title('Fourier Transform of reccovered object 90 deg', 'FontSize', fontSize, 'FontName', 'calibri')
colormap default,colorbar;
figure;
imagesc(log(abs(objectRecoverFTVV).*CTF_synth)); axis image
title('Fourier Transform of reccovered object 135 deg', 'FontSize', fontSize, 'FontName', 'calibri')
colormap default,colorbar;

pause
%% Recovery of delta and theta parameters

I0=abs(objectRecoverHH).^2;
I45=abs(objectRecoverHV).^2;
I90=abs(objectRecoverVH).^2;
I135=abs(objectRecoverVV).^2;

V1= (I90-I0)./(I90+I0);
V2= (I45-I135)./(I45+I135);

delta_recover = asin(sqrt(V1.^2+V2.^2));


theta_recover=zeros(m);
% for tt=1:m
%     for nn=1:m
%         if V2(tt,nn)>0
%             theta_recover(tt,nn) =  0.5*atan(V1(tt,nn)./V2(tt,nn));
%         else
%             theta_recover(tt,nn) =  0.5*atan(V1(tt,nn)./V2(tt,nn))+pi/2;
%         end
% 
%     end
% end

theta_recover = 0.5*asin(V1./sqrt(V1.^2+V2.^2));

figure;imagesc(abs(delta_recover));colorMap = jet(256);
colorMap = jet(256);  % Apply the colormap
colorbar;shading interp; 
axis image; title('recovered delta, rad')
caxis([0 3.14])
figure;imagesc(abs(theta_recover));colorMap = jet(256);colorbar; axis image;title('recovered theta,rad')
shading interp; caxis([0 3.14])
