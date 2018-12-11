clear
close all
tic
outerloop = 11;

loadfilename = sprintf('_11_oly_200.mat');
load(loadfilename)

x = downsample(downsample(x,4).',4).';
y = downsample(downsample(y,4).',4).';
exx = downsample(downsample(exx,4).',4).';
eyy = downsample(downsample(eyy,4).',4).';
u = downsample(downsample(u,4).',4).';
v = downsample(downsample(v,4).',4).';
% u(isnan(u)) = 0;
% v(isnan(v)) = 0;
% u = inpaintn(u,100);
% v = inpaintn(v,100);
% h1 = fspecial('average', 60);
% u=imfilter(u,h1);
% h1 = fspecial('average', 60);
% v=imfilter(v,h1);



DICx = x.*0.0827;
DICy = y.*0.0827;
%
% [newcoord1,newcoord2] = registerimages(DICx,DICy,exx,EBSDx,EBSDy,EBSDFID);
% save('registeredcoordinates1.mat','newcoord1')
% save('registeredcoordinates2.mat','newcoord2')
load('registeredcoordinates1.mat')
load('registeredcoordinates2.mat')

% [newcoord1_DIC,newcoord2_DIC] = registerimages_DIC(newcoord1.',newcoord2.',GrainID,DICx,DICy,DICexx);
% save('registeredcoordinates1_DIC.mat','newcoord1_DIC')
% save('registeredcoordinates2_DIC.mat','newcoord2_DIC')
load('registeredcoordinates1_DIC.mat')
load('registeredcoordinates2_DIC.mat')

GBdist = h5read('/Users/b119user/Downloads/Ajey/Al_Feb2018/20180109_ts3Al_07/GBNormals_v1.1/output/ANGfiles_Alroomtemp.h5','/GBdistance');
GrainID = h5read('/Users/b119user/Downloads/Ajey/Al_Feb2018/20180109_ts3Al_07/GBNormals_v1.1/output/ANGfiles_Alroomtemp.h5','/GID');
GBcoord = [newcoord1(GBdist==0),newcoord2(GBdist==0)];
GBID = GrainID(GBdist==0);

% quiver(GBcoord(:,1),GBcoord(:,2),GBnormals_GB(:,1),GBnormals_GB(:,2))

GrainIDGB = GrainID(GBdist==0);
GBidx = knnsearch(GBcoord,GBcoord,'K',50);
GrainIDGB_2 = GrainIDGB(GBidx);

% [p,q] = meshgrid(1:max(GrainID(:)),1:max(GrainID(:)));
% mask   = triu(ones(max(GrainID(:))), 1) > 0.5;
% pairs  = [p(mask) q(mask)];
countGBunique = 1;
for ii = 1:length(GBidx)
    GBpairID(ii,1) = GrainIDGB_2(ii,1);
    buff = find(GrainIDGB_2(ii,1:end)~=GrainIDGB_2(ii,1),1,'first');
    GBpairID(ii,2) = GrainIDGB_2(ii,buff);
    if (ii==1)
        GBID_pair(1,:) = GBpairID(ii,:);
        %     elseif (ismember(GBpairID(ii,:),GBID,'rows') || ismember(GBpairID(ii,:),fliplr(GBID),'rows'))
        %         continue
    else
        countGBunique = countGBunique + 1;
        GBID_pair(countGBunique,:) = GBpairID(ii,:);
    end
    if (ismember(GBpairID(ii,:),GBID_pair,'rows'))
        [~,Locb] = ismember(GBpairID(ii,:),GBID_pair,'rows');
        GBID_all(ii,1) = GBID_pair(Locb);
    elseif (ismember(GBpairID(ii,:),fliplr(GBID_pair),'rows'))
        [Lia,Locb] = ismember(GBpairID(ii,:),fliplr(GBID_pair),'rows');
        GBID_all(ii,1) = GBID_pair(Locb);
    end
    
end


GBunique = unique(GBID_pair,'rows');
%%
count = 0;
for ii = 1:length(GBunique)
    
    if (ii>1 && ismember(fliplr(GBunique(ii,:)),GBunique(1:ii-1,:),'rows'))
        continue
    else
%         figure('pos',[10,10,1000,800])
        count = count + 1;
        GBuniqueidx = union(find(all(bsxfun(@eq,GBID_pair,(GBunique(ii,:))), 2)),find(all(bsxfun(@eq,GBID_pair,fliplr(GBunique(ii,:))), 2)));
        GBuniquecoords{count} = GBcoord(GBuniqueidx,:);
    end
end
% 
% slidingbinary = zeros(42,1);
% % slidingbinary([18,19,20,25,33,34,35,36,40]) = 1;
% slidingbinary([1,11,18,19,20,25,27,33,34,36,40]) = 1;

slidingbinary = zeros(48,1);
slidingbinary([24,40,17,47,19]) = 1;
% exx(isnan(exx)) = 0;

%%
figure('pos',[10,10,2000,500])
surf(DICx,DICy,exx,'LineStyle','None')
view(2)
caxis([0,0.08])
colormap jet
axis equal
cb1 = colorbar;
hold on
for ii = 1:count
    if(slidingbinary(ii)==1)
        plot3(GBuniquecoords{ii}(:,1),GBuniquecoords{ii}(:,2),ones(size(GBuniquecoords{ii},1)).*1000,'k.','MarkerSize',20)
    end
end
xticks([])
yticks([])
% title('Sliding Boundaries')
set(gca,'FontSize',40)
xlim([min(DICx(:)),max(DICx(:))])
ylim([min(DICy(:)),max(DICy(:))])

title(cb1,'$\varepsilon_{xx}$','interpreter','latex','FontSize',50)
% saveas(gcf,'SlidingBoundaries.png')



toc
