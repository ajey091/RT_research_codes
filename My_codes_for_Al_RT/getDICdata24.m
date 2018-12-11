clear
close all

tic
loadfilename = sprintf('_11_oly_200.mat');
load(loadfilename)

exx(isnan(exx)) = 0;
eyy(isnan(eyy)) = 0;
exy(isnan(exy)) = 0;
exx = downsample(downsample(exx,4).',4).';
eyy = downsample(downsample(eyy,4).',4).';
exy = downsample(downsample(exy,4).',4).';
x = downsample(downsample(x,4).',4).';
y = downsample(downsample(y,4).',4).';

boxleftx = 339.5./0.0827;
boxrightx = 478.4./0.0827;
boxboty = 385.4./0.0827;
boxtopy = 420.1./0.0827;

figure('pos',[10,10,1200,1000])
surf(x.*0.0827,y.*0.0827,exx,'LineStyle','None')
view(2)
caxis([0,0.1])
xlim([boxleftx,boxrightx].*0.0827)
ylim([boxboty,boxtopy].*0.0827)
colorbar('eastoutside')
xticks([])
yticks([])
title('e_{xx} from DIC')
set(gca,'FontSize',28)
colormap jet
saveas(gcf,'GB24exx.png')

boxdim = mintersect(find(x>boxleftx),find(x<boxrightx),find(y>boxboty),find(y<boxtopy));

boxx = x(boxdim);
boxy = y(boxdim);

boxexx = mean(exx(boxdim))
boxeyy = mean(eyy(boxdim))
boxexy = mean(exy(boxdim))


% mean(boxexx(boxx==max(boxx)))
% mean(boxeyy(boxy==max(boxy)))
% mean(boxexy(boxy==max(boxy)))


toc