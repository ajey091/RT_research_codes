clear
close all force

tic

load('_11_oly_200.mat')
numbins = 50;
dsval = 4;
x = downsample(downsample(x,dsval).',dsval).';
y = downsample(downsample(y,dsval).',dsval).';
exx = downsample(downsample(exx,dsval).',dsval).';
exx = inpaintn(exx,100);
h1 = fspecial('average', 50);
exx=imfilter(exx,h1);
% exx = medfilt2(exx,[5,5]);
[~,edges] = histcounts(x(:),numbins);

ymid = mean(y(:));

f = waitbar(0,'Running');
for ii = 1:numbins
    waitbar(ii/numbins)
    coordofinterest(ii,:) = [(edges(ii+1)+edges(ii))/2,ymid];
    idx = knnsearch([x(:),y(:)],coordofinterest(ii,:));
    nearestmeshpoint(ii,:) = [x(idx),y(idx)];
    strain11lineDIC(ii,:) = exx(idx);
%     nonlocalmeshpoints =  rangesearch([x(:),y(:)],nearestmeshpoint(ii,:),40);
%     nonlocalDICstrain(ii,:) = mean(exx(nonlocalmeshpoints{:}));
end
csvfilename = sprintf('/Users/b119user/Downloads/Ajey/Al_Ajey/Aug2017/Test_specimen/ANGfiles/PaddedANGfiles/STLfiles/Results_extrudedpad/Results_4_16/Al_full_model_lowres_11_converged_3/Al_full_model_11_all.csv');
% csvfilename = sprintf('/Users/b119user/Downloads/Ajey/Al_Feb2018/20180109_ts3Al_07/ANGfiles/Padded/STLfiles/Results/Al_full_model_lowres_highT_07_mantle_size15/Al_lowres_model_mantle15.csv');

data = csvread(csvfilename,1,0);
nodescoor = data(:,8:10);
stress11 = data(:,1);

% [newcoord1,newcoord2] = registerimages_DIC(x,y,exx,nodescoor(:,1),nodescoor(:,2),stress11);
% save('registeredcoordinates1_DICCP.mat','newcoord1')
% save('registeredcoordinates2_DICCP.mat','newcoord2')
load('registeredcoordinates1_DICCP.mat')
load('registeredcoordinates2_DICCP.mat')

nodescoor2 = [newcoord1,newcoord2];
figure,
plot(nodescoor(:,1),nodescoor(:,2),'.k')
hold on
stress22 = data(:,2);
stress33 = data(:,3);
stress12 = data(:,4);
stress13 = data(:,5);
stress23 = data(:,6);
strain11 = data(:,7);





% [~,edges] = histcounts(nodescoor2(:,1),numbins);
% ymid = mean(nodescoor2(:,2));
for ii = 1:numbins
    
    %     coordofinterest(ii,:) = [(edges(ii+1)+edges(ii))/2,ymid];
    %     idx = knnsearch(nodescoor2(:,1:2),coordofinterest(ii,:));
    idx = knnsearch(nodescoor2,coordofinterest(ii,:));
%     nearestmeshpoint(ii,:) = nodescoor2(idx,1:2);
%     nonlocalmeshpoints =  rangesearch(nodescoor2(:,1:2),nearestmeshpoint(ii,:),0);
%     nonlocalavgstress11(ii,:) = mean(stress11(nonlocalmeshpoints{:}));
%     stress11line(ii) = stress11(idx);
    strain11line(ii) = strain11(idx);
%     nonlocalcoords{ii} = nodescoor2(nonlocalmeshpoints{:},:);
    
end

%%
figure('pos',[10,10,1000,600]),
plot(linspace(0,1,numbins-2),strain11line(2:end-1),'o-','LineWidth',4,'Color',[0.4,0,1]);
hold on;
plot(linspace(0,1,numbins-2),strain11lineDIC(2:end-1),'o-','LineWidth',4,'Color',[0.8,0.2,0])
% plot(linspace(0,1,numbins-2),strain11line15(2:end-1),'o-','LineWidth',4)
legend('CPFE','DIC','location','best')
% ylim([-0.01,0.15])
% title('Strain in the loading direction')
xlabel('Normalized distance along specimen length')
ylabel('Strain in the loading direction')
set(gca,'YAxisLocation','right')
set(gca,'FontSize',36)
saveas(gcf,'DICCPFEStraincompare.png')


toc