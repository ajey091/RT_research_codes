clear
close all

load('_11_oly_200.mat')
EBSDdata = load('EBSD_euler_aligned_with_sample.mat');
exx = downsample(downsample(exx,3).',3).';
x = downsample(downsample(x,3).',3).';
y = downsample(downsample(y,3).',3).';
figure('pos',[100 100 1800 800])
surf(x,y,exx,'LineStyle','None');view(2)
hold on

rectangle('Position',[x(180,1960) y(180,1960) x(520,3080)-x(180,1960) y(520,3080)-y(180,1960)],...
    'EdgeColor','r',  'LineWidth',3)
text(x(180,1960)*0.9+(x(520,3080)-x(180,1960))/2,y(180,1960)*0.2+(y(520,3080)),'Region 1','FontSize',16,'color','w')

rectangle('Position',[x(405,60) y(405,60) x(580,550)-x(405,60) y(580,550)-y(405,60)],...
    'EdgeColor','r',  'LineWidth',3)
text(x(405,60)*2,y(405,60)*0.08+(y(580,550)),'Region 2','FontSize',16,'color','w')

rectangle('Position',[x(78,25) y(78,25) x(265,525)-x(78,25) y(265,525)-y(78,25)],...
    'EdgeColor','r',  'LineWidth',3)
text(x(78,25)*5,y(78,25)*0.3+(y(265,525)),'Region 3','FontSize',16,'color','w')

rectangle('Position',[x(500,3550) y(500,3550) x(780,4580)-x(500,3550) y(780,4580)-y(500,3550)],...
    'EdgeColor','r',  'LineWidth',3)
text(x(500,3550)*0.97+(x(780,4580)-x(500,3550))/2,y(500,3550)*0.05+(y(780,4580)),'Region 4','FontSize',16,'color','w')

rectangle('Position',[x(620,4640) y(620,4640) x(760,4950)-x(620,4640) y(760,4950)-y(620,4640)],...
    'EdgeColor','r',  'LineWidth',3)
text(x(620,4640)*0.963+(x(760,4950)-x(620,4640))/2,y(620,4640)*0.05+(y(760,4950)),'Region 5','FontSize',16,'color','w')

xlim([min(x(:)),max(x(:))])
caxis([0 0.06])
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gcf,'color','w')
colorbar
% export_fig('-png','-nocrop','-opengl','-r200','Regions_all')

slipplanes=1/sqrt(3)*[1 1 -1; % Slip planes
                      1 -1 -1; 
                      1 -1  1;
                      1  1  1];
        

Region1_exx = exx(180:520,1960:3080);
Region1_exx = inpaintn(Region1_exx,100);
figure('pos',[100 100 1800 800])
surf(Region1_exx,'LineStyle','None');view(2)
caxis([0 0.06])
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gcf,'color','w')
title('Region 1','FontSize',20)
% export_fig('-png','-nocrop','-opengl','-r200','Region_1')

Region2_exx = exx(405:580,60:550);
Region2_exx = inpaintn(Region2_exx,100);
figure('pos',[100 100 1800 800])
surf(Region2_exx,'LineStyle','None');view(2)
caxis([0 0.06])
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gcf,'color','w')
title('Region 2','FontSize',20)
% export_fig('-png','-nocrop','-opengl','-r200','Region_2')

Region3_exx = exx(78:265,25:525);
Region3_exx = inpaintn(Region3_exx,100);
figure('pos',[100 100 1800 800])
surf(Region3_exx,'LineStyle','None');view(2)
caxis([0 0.06])
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gcf,'color','w')
title('Region 3','FontSize',20)
% export_fig('-png','-nocrop','-opengl','-r200','Region_3')

Region4_exx = exx(500:780,3550:4580);
Region4_exx = inpaintn(Region4_exx,100);
figure('pos',[100 100 1800 800])
surf(Region4_exx,'LineStyle','None');view(2)
caxis([0 0.06])
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gcf,'color','w')
title('Region 4','FontSize',20)
% export_fig('-png','-nocrop','-opengl','-r200','Region_4')

Region5_exx = exx(620:760,4640:4950);
Region5_exx = inpaintn(Region5_exx,100);
figure('pos',[100 100 1800 800])
surf(Region5_exx,'LineStyle','None');view(2)
caxis([0 0.06])
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gcf,'color','w')
title('Region 5','FontSize',20)
% export_fig('-png','-nocrop','-opengl','-r200','Region_5')
