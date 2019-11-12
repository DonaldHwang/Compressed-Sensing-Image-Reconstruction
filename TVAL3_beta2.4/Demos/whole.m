function whole

pa = genpath(pwd);% 获得文件夹data下所有子文件的路径，这些路径存在字符串p中，以';'分割
length_p = size(pa,2);%字符串p的长度
path = {};%建立一个单元数组，数组的每个单元中包含一个目录
temp = [];
for i = 1:length_p %寻找分割符';'，一旦找到，则将路径temp写入path数组中
    if pa(i) ~= ';'
        temp = [temp pa(i)];
    else 
        temp = [temp '\']; %在路径的最后加入 '\'
        path = [path ; temp];
        temp = [];
    end
end  
clear pa length_p temp;
%至此获得data文件夹及其所有子文件夹（及子文件夹的子文件夹）的路径，存于数组path中。
%下面是逐一文件夹中读取图像
file_num = size(path,1);% 子文件夹的个数
for i = 1:file_num
    file_path =  path{i}; % 图像文件夹路径
    dest_path = strcat(path{i});
    img_path_list = dir(strcat(file_path,'*.jpg'));
    img_num = length(img_path_list); %该文件夹中图像数量
    if img_num > 0
        for j = 1:img_num
            image_name = img_path_list(j).name;% 图像名
            rec_image_name = strcat(erase(image_name,".jpg"),"_rec.png");
            image =  imread(strcat(file_path,image_name));
            fprintf('%d %d %s\n',i,j,strcat(file_path,image_name));% 显示正在处理的路径和图像名
            result = Copy_of_demo_barbara(image);
            imwrite(result,strcat(file_path,rec_image_name));
        end
    end
end


% clear all; close all;
% 
% fullscreen = get(0,'ScreenSize');
% 
% % problem size
% ratio = .25;
% sidelength = 32;
% N = sidelength^2;
% M = round(ratio*N);
% 
% % original image
% bbr = importdata('0_30.jpg');
% Im1 = double(bbr(:,:,1));
% Im2 = double(bbr(:,:,2));
% Im3 = double(bbr(:,:,3));
% 
% % generate measurement matrix
% p = randperm(N);
% picks = p(1:M);
% for ii = 1:M
%     if picks(ii) == 1
%         picks(ii) = p(M+1);
%         break;
%     end
% end
% perm = randperm(N); % column permutations allowable
% A = @(x,mode) dfA(x,picks,perm,mode);
% 
% % observation
% b1 = A(Im1(:),1);
% bavg1 = mean(abs(b1));
% b2 = A(Im2(:),1);
% bavg2 = mean(abs(b2));
% b3 = A(Im3(:),1);
% bavg3 = mean(abs(b3));
% 
% % add noise
% sigma = 0.04;  % noise std
% noise = randn(M,1);
% b1 = b1 + sigma*bavg1*noise;
% b2 = b2 + sigma*bavg2*noise;
% b3 = b3 + sigma*bavg3*noise;
% 
% % set the optional paramaters
% clear opts
% opts.mu = 2^12;
% opts.beta = 2^6;
% opts.mu0 = 2^4;       % trigger continuation shceme
% opts.beta0 = 2^-2;    % trigger continuation shceme
% opts.maxcnt = 10;
% opts.tol_inn = 1e-3;
% opts.tol = 1E-6;
% opts.maxit = 300;
% 
% % reconstruction
% t = cputime;
% [estIm1, out1] = TVAL3(A,b1,sidelength,sidelength,opts);
% estIm1 = estIm1 - min(estIm1(:));
% [estIm2, out2] = TVAL3(A,b2,sidelength,sidelength,opts);
% estIm2 = estIm2 - min(estIm2(:));
% [estIm3, out3] = TVAL3(A,b3,sidelength,sidelength,opts);
% estIm3 = estIm3 - min(estIm3(:));
% 
% 
% t = cputime - t;
% re_er = norm(estIm1-Im1,'fro')/norm(Im1,'fro');
% 
% s1 = uint8(estIm1);
% s2 = uint8(estIm2);
% s3 = uint8(estIm3);
% figure('Name','Barbara','Position',...
%     [fullscreen(1) fullscreen(2) fullscreen(3) fullscreen(4)]);
% colormap(gray);
% estIm = cat(3,s1,s2,s3);
% subplot(1,2,1);
% imshow(bbr,[]);
% title(sprintf('Original %dx%d Barbara',sidelength,sidelength),'fontsize',16);
% subplot(1,2,2);
% imshow(estIm,[]);
% title(sprintf('Recovered Barbara with %2.0f%% measurements',ratio*100),'fontsize',16);
% xlabel(sprintf('Noise level: %2d%%  \n Rel-Err: %4.2f%%,   CPU time: %4.2fs',100*sigma,100*re_er,t),'fontsize',14);
% 
% 
% plotting = 0;
% if plotting
%     figure(2);
%     subplot(241); plot(out.lam1); title('\_al: ||w||');
%     subplot(242); plot(out.lam2); title('\_al: ||Du-w||^2');
%     subplot(243); plot(out.lam3); title('\_al: ||Au-f||^2');
%     subplot(244); plot(abs(out.obj),'b-'); title('\_al: objective values');
%     subplot(245); plot(out.res); title('\_al: residue');
%     subplot(246); plot(abs(out.tau)); title('\_al: steplenths');
%     subplot(247); plot(out.itrs); title('\_al: inner iterations');
%     subplot(248); plot(abs(out.C),'r-'); title('\_al: reference vlaues');
%     
%     figure(3);
%         semilogy(1:length(out.lam1),out.lam1,'b*:',1:length(out.lam2),sqrt(out.lam2),'rx:',...
%         1:length(out.lam3),sqrt(out.lam3),'g.--', 1:length(out.f),sqrt(out.f),'m+-');
%     legend('lam1(||w||_1)','lam2(||D(d_tu)-w||_2)','lam3(||Au-b||_2)','obj function');
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % dfA
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y = dfA(x,picks,perm,mode)
% switch mode
%     case 1
%         y = A_fWH(x,picks,perm);
%     case 2
%         y = At_fWH(x,picks,perm);
%     otherwise
%         error('Unknown mode passed to f_handleA!');
% end
