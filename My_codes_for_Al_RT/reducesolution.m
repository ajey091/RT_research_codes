clear
close all

tic
for ii = 1:11
    loadfilename = sprintf('_%.2d_oly_200.mat',ii);
    load(loadfilename)
    
    % exx_ds = downsample(downsample(exx,4).',4).';
    exx(isnan(exx)) = 0;
    % exx_ds = medfilt2(exx_ds);
    h1 = fspecial('average', 250);
    buff2=imfilter(exx,h1);
    
    %
    % figure('pos',[10,10,1800,600])
    % surf(exx_ds,'LineStyle','None')
    % view(2)
    % caxis([0,0.15])
    
    figure('pos',[10,10,2200,600])
    surf(x,y,buff2,'LineStyle','None')
    view(2)
    caxis([0,0.04])
    xlim([min(x(:)),max(x(:))])
    colorbar('southoutside')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title('e_{xx} from DIC')
    set(gca,'FontSize',28)
    colormap jet
    savefilename = sprintf('lowresDIC_%.2d.png',ii);
    saveas(gcf,savefilename)
    close all
end

toc