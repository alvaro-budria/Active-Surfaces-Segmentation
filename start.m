%close all;
clearvars;
clc

% %%% circles grayscale%%%
[I, ni, nj] = read_image('circles.png', true, false, true, false);
mu=1; nu=0.5;
lambda1=1; lambda2=1;
epHeaviside=1; eta=0.01;
tol=0.000001;
dt=5*(10^0)/mu; iterMax=350;
reIni=100000;
init="cone";
directory="circles_init";
% iters: 14
% default:
% commments:
%%%%%%%%%%%%%%%%%%%%


% %%% Blobs grayscale%%%
% [I, ni, nj] = read_image('blobs.png', true, false, false, false);
% I = (1-(1 - I | binornd(1, 0.2, ni, nj)) | binornd(1, 0.2, ni, nj));
% mu=1; nu=-0.0001;
% lambda1=10; lambda2=10;
% epHeaviside=1; eta=0.01;
% tol=0.000001;
% dt=5*(10^0)/mu; iterMax=350;
% reIni=100000;
% init="checkerboard";
% directory="blobs_noise";
% % iters: 14
% % default:
% % commments:
% %%%%%%%%%%%%%%%%%%%%%


% %%% Noised Circles %%%
% [I, ni, nj] = read_image('./images/noisedCircles.tif', true, false, true, false);
% mu=0.125; nu=0.085;
% %mu=0.1;
% lambda1=0.01; lambda2=0.75;
% epHeaviside=1; eta=0.01;
% tol=0.05;
% dt=5*(10^7)/mu; iterMax=350;
% reIni=100000;
% init="cone";
% directory="noisedCircles";
% % iters: 14
% % default:
% % commments:
% %%%%%%%%%%%%%%%%%%%%%

% %%% Phantom 17 %%%
% [I, ni, nj] = read_image('./images/phantom17.bmp', true, false, true, true);
% mu=1; nu=0.1;
% %mu=0.1;
% lambda1=2; lambda2=2;
% epHeaviside=1; eta=0.01;
% tol=0.05;
% dt=5*(10^-1)/mu; iterMax=350;
% reIni=100000;
% init="cone";
% directory="phantom17";
% % iters: 14
% % default:
% % commments:
% %%%%%%%%%%%%%%%%%%%%%

% %%% Europe Binary %%%
% [I, ni, nj] = read_image('./images/map_europe_binary.png', false, false, true, true);
% mu=0.005; nu=0.01;
% lambda1=1; lambda2=1;
% epHeaviside=1; eta=0.01;
% tol=0.1;
% dt=5*(10^-1)/mu; iterMax=30;
% reIni=15;
% init="checkerboard";
% directory=strcat("europe_binary", "_", num2str(mu), "_", num2str(nu), "_", num2str(lambda1), "_", num2str(lambda2), "_", num2str(init), "_", num2str(dt));
% % iters: 14
% % default:
% % commments:
% %%%%%%%%%%%%%%%%%%%%%

% %%% Europe Grayscale %%%
% [I, ni, nj] = read_image('./images/map_europe.png', true, false, true, true);
% mu=0.001; nu=0.0;
% lambda1=1; lambda2=0.5;
% epHeaviside=1; eta=0.01;
% tol=0.1;
% dt=5*(10^-1)/mu; iterMax=100;
% reIni=20;
% init="checkerboard";
% directory="europe_gray";
% % iters: 14
% % default:
% % commments:
% %%%%%%%%%%%%%%%%%%%%

% %%% World Binary %%%
% [I, ni, nj] = read_image('./images/map_world_binary.png', false, false, true, true);
% mu=0.005; nu=0.035;
% lambda1=1; lambda2=1;
% epHeaviside=1; eta=0.01;
% tol=0.1;
% dt=5*(10^-10)/mu; iterMax=30;
% reIni=3;
% init="custom_world";
% directory="world_binary_reinit";
% % iters: 14
% % default:
% % commments:
% %%%%%%%%%%%%%%%%%%%%%

% %%% Face1 Grayscale %%%
% [I, ni, nj] = read_image('./images/face1.png', true, false, true, true);
% mu=1.75; nu=-0.5;
% lambda1=1; lambda2=1;
% epHeaviside=1.5; eta=0.01;
% tol=0.1;
% dt=5*(10^-0)/mu; iterMax=245;
% reIni=20;
% init="checkerboard";
% directory="face1";
% % iters: 14
% % default:
% % commments:
% %%%%%%%%%%%%%%%%%%%%


% %%% Face6 Grayscale %%%
% [I, ni, nj] = read_image('./images/face6.png', true, false, true, true);
% mu=5.3; nu=-1.2;
% lambda1=4; lambda2=0.001;
% epHeaviside=1; eta=0.01;
% tol=0.1;
% dt=5*(10^1)/mu; iterMax=100;
% reIni=45;
% init="checkerboard";
% directory="face6";
% % iters: 14
% % default:
% % commments:
% %%%%%%%%%%%%%%%%%%%%


%%% Face6 Color %%%
% isColor = true;
% [I, ni, nj] = read_image('./images/face6.png', true, isColor, true, true);
% mu=1; nu=-3.0;
% lambda1=40; lambda2=0.01;
% epHeaviside=1; eta=0.01;
% tol=0.1;
% dt=5*(10^1)/mu; iterMax=100;
% reIni=45;
% init="checkerboard";
% directory=strcat("face6", "_isColor-", string(isColor));
% % iters: 14
% % default:
% % commments:
% %%%%%%%%%%%%%%%%%%%%


% %%% Landscape1 Grayscale %%%
% [I, ni, nj] = read_image('./images/landscape1.png', true, false, true, true);
% mu=13; nu=-0.5;
% lambda1=10; lambda2=0.01;
% epHeaviside=1; eta=0.01;
% tol=0.0000001;
% dt=1*(10^1)/mu; iterMax=100;
% reIni=40;
% init="checkerboard";
% directory="landscape1";
% % iters: 14
% % default:
% % commments:
% %%%%%%%%%%%%%%%%%%%%


% %%% Landscape2 Grayscale %%%
% [I, ni, nj] = read_image('./images/landscape3.png', true, false, true, true);
% mu=100; nu=-0.5;
% lambda1=20; lambda2=0.01;
% epHeaviside=1.5; eta=0.01;
% tol=0.001;
% dt=5*(10^-1)/mu; iterMax=120;
% reIni=40;
% init="checkerboard";
% directory="landscape3";
% % iters: 14
% % default:
% % commments:
% %%%%%%%%%%%%%%%%%%%%


% %%% 00015 Grayscale %%%
% [I, ni, nj] = read_image('./images/00015.jpg', true, false, true, true);
% mu=10; nu=-1;
% lambda1=10; lambda2=0.1;
% epHeaviside=11; eta=0.01;
% tol=0.001;
% dt=1*(10^-0)/mu; iterMax=120;
% reIni=40;
% init="checkerboard";
% directory="00015";
% % iters: 14
% % default:
% % commments:
% %%%%%%%%%%%%%%%%%%%%


% %%% 00016 Grayscale %%%
% [I, ni, nj] = read_image('./images/00016.jpg', true, false, true, true);
% mu=0.5; nu=-0.3;
% lambda1=10; lambda2=0.1;
% epHeaviside=11; eta=0.01;
% tol=0.001;
% dt=5*(10^0)/mu; iterMax=120;
% reIni=40;
% init="checkerboard";
% directory="00016";
% % iters: 14
% % default:
% % commments:
% %%%%%%%%%%%%%%%%%%%%

% %%% 00022 Grayscale %%%
% [I, ni, nj] = read_image('./images/00022.jpg', true, false, true, true);
% mu=0.1; nu=-0.4;
% lambda1=10; lambda2=0.1;
% epHeaviside=11; eta=0.01;
% tol=0.001;
% dt=5*(10^0)/mu; iterMax=120;
% reIni=40;
% init="checkerboard";
% directory="0002";
% % iters: 14
% % default:
% % commments:
% %%%%%%%%%%%%%%%%%%%%


% %%% phantom 18 Grayscale %%%
% [I, ni, nj] = read_image('./images/phantom18.bmp', true, false, true, false);
% mu=0.00000001; nu=-0.01;
% lambda1=0.001; lambda2=0.001;
% epHeaviside=0.1; eta=0.1;
% tol=0.001;
% dt=5*(10^-9)/mu; iterMax=120;
% reIni=40;
% init="cone";
% directory="phantom18";
% % iters: 14
% % default:
% % commments:
% %%%%%%%%%%%%%%%%%%%%


% %% eagle Grayscale %%%
% [I, ni, nj] = read_image('./images/eagle.png', true, false, true, false);
% mu=0.3; nu=-0.12;
% lambda1=100; lambda2=0.1;
% epHeaviside=3; eta=0.1;
% tol=0.001;
% dt=5*(10^-1)/mu; iterMax=120;
% reIni=40;
% init="cone";
% directory="eagle_gray";
% iters: 14
% default:
% commments:
% %%%%%%%%%%%%%%%%%%%

% 
% %%% reja Grayscale %%%
% [I, ni, nj] = read_image('./images/reja.png', true, false, true, false);
% mu=0.17; nu=-2.5;
% lambda1=15; lambda2=0.1;
% epHeaviside=3; eta=0.1;
% tol=0.001;
% dt=5*(10^-1)/mu; iterMax=120;
% reIni=40;
% init="cone";
% directory="reja_gray";
% % iters: 14
% % default:
% % commments:
% %%%%%%%%%%%%%%%%%%%%


%lambda1=10^-3; %Hola carola problem
%lambda2=10^-3; %Hola carola problem

%eta=0.01;

%dt=(10^-2)/mu; 

%reIni=0; %Try both of them
%reIni=500;


directory=strcat(num2str(directory), "_", num2str(mu), "_", num2str(nu), "_", num2str(lambda1), "_", num2str(lambda2), "_", num2str(init), "_", num2str(dt));


[X, Y]=meshgrid(1:nj, 1:ni);

%%Initial phi
if init == "cone"
    phi_0 = (-sqrt( ( X-round(ni/2)).^2 + (Y-round(nj/2)).^2)+50);
elseif init == "checkerboard"
    phi_0 = sin(pi/5 * X) .* sin(pi/5 * Y);
elseif init == "custom_world"
    phi_0 = (-sqrt( ( X+round(ni/1)).^2 + (Y-round(nj/200)).^2)+100);
end

%%% This initialization allows a faster convergence for phantom 18
%phi_0=(-sqrt( ( X-round(ni/2)).^2 + (Y-round(nj/4)).^2)+50);
%Normalization of the initial phi to [-1 1]
%phi_0=phi_0-min(phi_0(:));
%phi_0=2*phi_0/max(phi_0(:));
%phi_0=phi_0-1;

%phi_0=I; %For the Hola carola problem

phi_0=phi_0-min(phi_0(:));
phi_0=2*phi_0/max(phi_0(:));
phi_0=phi_0-1;


%%Explicit Gradient Descent
seg=sol_ChanVeseIpol_GDExp( I, phi_0, mu, nu, eta, lambda1, lambda2, tol, epHeaviside, dt, iterMax, reIni, directory );
% imshow("setUnits", seg)

writerObj = VideoWriter(strcat(directory, '/chan_vese.mp4')); % new video
writerObj.FrameRate = 2;
open(writerObj);
imgs_names = dir(fullfile(directory, '/*.png'));
for K = 1 : numel(imgs_names)
    baseFileName = imgs_names(K).name;
    fullFileName = fullfile(imgs_names(K).folder, baseFileName);
    % filename = sprintf('a%04d.tif', K);
    thisimage = imread(fullFileName);
    writeVideo(writerObj, thisimage);
end
close(writerObj);
