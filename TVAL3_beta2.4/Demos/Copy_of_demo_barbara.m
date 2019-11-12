function estIm=Copy_of_demo_barbara(bbr)
%
% This demo shows how TVAL3 handles a complicated image and how to trigger
% the continuation scheme in order to speed up the convergence.
%
% I: 256x256 Barbara (real, two-dimentional)
% A: permuted Walsh Hadamard transform (real)
% f: observation with noise (real)
% 
% Written by: Chengbo Li
% Advisor: Prof. Yin Zhang and Wotao Yin 
% CAAM department, Rice University
% 05/21/2009


path(path,genpath(pwd));
% fullscreen = get(0,'ScreenSize');

% problem size
ratio = .25;
sidelength = 32;
N = sidelength^2;
M = round(11);

% original image

Im1 = double(bbr(:,:,1));
Im2 = double(bbr(:,:,2));
Im3 = double(bbr(:,:,3));

% generate measurement matrix
p = randperm(N);
picks = p(1:M);
for ii = 1:M
    if picks(ii) == 1
        picks(ii) = p(M+1);
        break;
    end
end
perm = randperm(N); % column permutations allowable
A = @(x,mode) dfA(x,picks,perm,mode);

% observation
b1 = A(Im1(:),1);
bavg1 = mean(abs(b1));
b2 = A(Im2(:),1);
bavg2 = mean(abs(b2));
b3 = A(Im3(:),1);
bavg3 = mean(abs(b3));

% add noise
sigma = 0.04;  % noise std
noise = randn(M,1);
b1 = b1 + sigma*bavg1*noise;
b2 = b2 + sigma*bavg2*noise;
b3 = b3 + sigma*bavg3*noise;

% set the optional paramaters
clear opts
opts.mu = 2^12;
opts.beta = 2^6;
opts.mu0 = 2^4;       % trigger continuation shceme
opts.beta0 = 2^-2;    % trigger continuation shceme
opts.maxcnt = 10;
opts.tol_inn = 1e-3;
opts.tol = 1E-6;
opts.maxit = 300;

% reconstruction
% t = cputime;
[estIm1, out1] = TVAL3(A,b1,sidelength,sidelength,opts);
estIm1 = estIm1 - min(estIm1(:));
[estIm2, out2] = TVAL3(A,b2,sidelength,sidelength,opts);
estIm2 = estIm2 - min(estIm2(:));
[estIm3, out3] = TVAL3(A,b3,sidelength,sidelength,opts);
estIm3 = estIm3 - min(estIm3(:));


% t = cputime - t;
% re_er = norm(estIm1-Im1,'fro')/norm(Im1,'fro');

s1 = uint8(estIm1);
s2 = uint8(estIm2);
s3 = uint8(estIm3);
% figure('Name','Barbara','Position',...
%     [fullscreen(1) fullscreen(2) fullscreen(3) fullscreen(4)]);
% colormap(gray);
estIm = cat(3,s1,s2,s3);
% subplot(1,2,1);
% imshow(bbr,[]);
% title(sprintf('Original %dx%d Barbara',sidelength,sidelength),'fontsize',16);
% subplot(1,2,2);
% imshow(estIm,[]);
% title(sprintf('Recovered Barbara with %2.0f%% measurements',ratio*100),'fontsize',16);
% xlabel(sprintf('Noise level: %2d%%  \n Rel-Err: %4.2f%%,   CPU time: %4.2fs',100*sigma,100*re_er,t),'fontsize',14);


plotting = 0;
if plotting
    figure(2);
    subplot(241); plot(out.lam1); title('\_al: ||w||');
    subplot(242); plot(out.lam2); title('\_al: ||Du-w||^2');
    subplot(243); plot(out.lam3); title('\_al: ||Au-f||^2');
    subplot(244); plot(abs(out.obj),'b-'); title('\_al: objective values');
    subplot(245); plot(out.res); title('\_al: residue');
    subplot(246); plot(abs(out.tau)); title('\_al: steplenths');
    subplot(247); plot(out.itrs); title('\_al: inner iterations');
    subplot(248); plot(abs(out.C),'r-'); title('\_al: reference vlaues');
    
    figure(3);
        semilogy(1:length(out.lam1),out.lam1,'b*:',1:length(out.lam2),sqrt(out.lam2),'rx:',...
        1:length(out.lam3),sqrt(out.lam3),'g.--', 1:length(out.f),sqrt(out.f),'m+-');
    legend('lam1(||w||_1)','lam2(||D(d_tu)-w||_2)','lam3(||Au-b||_2)','obj function');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% dfA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = dfA(x,picks,perm,mode)
switch mode
    case 1
        y = A_fWH(x,picks,perm);
    case 2
        y = At_fWH(x,picks,perm);
    otherwise
        error('Unknown mode passed to f_handleA!');
end
