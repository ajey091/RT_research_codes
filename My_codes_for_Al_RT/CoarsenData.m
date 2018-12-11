%{
clear
load('_11_oly_200.mat')

%%
close all


u_ds = downsample(downsample(u,4).',4).';
u_ds(isnan(u_ds)) = 0;
u_ds = medfilt2(u_ds);
h1 = fspecial('average', 10);
buff2=imfilter(u_ds,h1);


figure
surf(u_ds,'LineStyle','None')
view(2)
caxis([-1000,1000])

figure('pos',[10,10,1800,600])
surf(buff2,'LineStyle','None')
view(2)
% caxis([-1000,1000])
colorbar
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'FontSize',28)

% saveas(gcf,'03disp.png')

clear
load('_06_oly_200.mat')

%%
close all


u_ds = downsample(downsample(u,4).',4).';
u_ds(isnan(u_ds)) = 0;
u_ds = medfilt2(u_ds);
h1 = fspecial('average', 10);
buff2=imfilter(u_ds,h1);


figure
surf(u_ds,'LineStyle','None')
view(2)
caxis([-1000,1000])

figure('pos',[10,10,1800,600])
surf(buff2,'LineStyle','None')
view(2)
% caxis([-1000,1000])
colorbar
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'FontSize',28)

% saveas(gcf,'06disp.png')


clear
% load('_09_oly_200.mat')

%%
close all


u_ds = downsample(downsample(u,4).',4).';
u_ds(isnan(u_ds)) = 0;
u_ds = medfilt2(u_ds);
h1 = fspecial('average', 10);
buff2=imfilter(u_ds,h1);


figure
surf(u_ds,'LineStyle','None')
view(2)
caxis([-1000,1000])

figure('pos',[10,10,1800,600])
surf(buff2,'LineStyle','None')
view(2)
% caxis([-1000,1000])
colorbar
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'FontSize',28)

% saveas(gcf,'09disp.png')
%}

clear
load('_11_oly_200.mat')

%%
close all

x_ds = downsample(downsample(x,4).',4).';
y_ds = downsample(downsample(y,4).',4).';

exx_ds = downsample(downsample(exx,4).',4).';
exx_ds(isnan(exx_ds)) = 0;
exx_ds = medfilt2(exx_ds);
h1 = fspecial('average', 50);
buff2=imfilter(exx_ds,h1);


% figure('pos',[10,10,1600,300])
% surf(x_ds,y_ds,exx_ds,'LineStyle','None')
% view(2)
% caxis([0,0.04])
% axis equal

%%
figure('pos',[10,10,1600,300])
surf(x_ds,y_ds,buff2,'LineStyle','None')
view(2)
caxis([0,0.04])
axis equal
h = colorbar('westoutside');
title(h,'$\varepsilon_{xx}$ (DIC)','Interpreter','Latex','FontSize',50)
colormap jet
xlim([min(x_ds(:)),max(x_ds(:))])
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'FontSize',34)
set(gca,'XColor','none')
set(gca,'YColor','none')
% saveas(gcf,'11disp.png')
