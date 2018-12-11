clear
close all
tic
for outerloop = 11
   
    loadfilename = sprintf('_%2.2d_oly_200.mat',outerloop);
    load(loadfilename)
    
    x = downsample(downsample(x,4).',4).';
    y = downsample(downsample(y,4).',4).';
    exx = downsample(downsample(exx,4).',4).';
    eyy = downsample(downsample(eyy,4).',4).';
    u = downsample(downsample(u,4).',4).';
    v = downsample(downsample(v,4).',4).';
    % u(isnan(u)) = 0;
    % v(isnan(v)) = 0;
    u = inpaintn(u,100);
    v = inpaintn(v,100);
    % h1 = fspecial('average', 60);
    % u=imfilter(u,h1);
    % h1 = fspecial('average', 60);
    % v=imfilter(v,h1);
    for ii = 1:length(u)
        u_corrected (:,ii) = u(:,ii) - mean(u(:,ii));
        
    end
    for ii = 1:size(v,1)
        
        v_corrected (ii,:) = v(ii,:) - mean(v(ii,:));
    end
    
    
    DICx = x.*0.0827;
    DICy = y.*0.0827;
    EBSDx = squeeze(h5read('Test_specimen_output.dream3d','/DataContainers/ImageDataContainer/CellData/X Position'));
    EBSDy = squeeze(h5read('Test_specimen_output.dream3d','/DataContainers/ImageDataContainer/CellData/Y Position'));
    EBSDFID = squeeze(h5read('Test_specimen_output.dream3d','/DataContainers/ImageDataContainer/CellData/FeatureIds'));
    EBSDFID = flipud(EBSDFID);
    %
%     [newcoord1,newcoord2] = registerimages(DICx,DICy,exx,EBSDx,EBSDy,EBSDFID);
%     save('registeredcoordinates1.mat','newcoord1')
%     save('registeredcoordinates2.mat','newcoord2')
    loadfilename1 = sprintf('registeredcoordinates1_%d.mat',outerloop);
    loadfilename2 = sprintf('registeredcoordinates2_%d.mat',outerloop);
    load(loadfilename1)
    load(loadfilename2)
    
    GBdist = h5read('/Users/b119user/Downloads/Ajey/Al_Feb2018/20180109_ts3Al_07/GBNormals_v1.1/output/ANGfiles_Alroomtemp.h5','/GBdistance');
    GBdist(GBdist~=0) = NaN;
    GBdist(GBdist==0) = 100;
    % figure('pos',[10,10,1800,600])
    % s1 = surf(DICx,DICy,exx,'LineStyle','None');
    % hold on
    % s2 = surf(newcoord1,newcoord2,GBdist,'LineStyle','None');
    % view(2)
    % s1.FaceAlpha = 0.5;
    % s2.FaceAlpha = 0.5;
    
    
    GBnormals = h5read('/Users/b119user/Downloads/Ajey/Al_Feb2018/20180109_ts3Al_07/GBNormals_v1.1/output/ANGfiles_Alroomtemp.h5','/GBNormals');
    GBnormals = permute(GBnormals,[2,3,1]);
    GBnormals = reshape(GBnormals,[size(EBSDx,1)*size(EBSDx,2),3]);
    GBTangents = cross(GBnormals,repmat([0,0,1],[size(EBSDx,1)*size(EBSDx,2),1]),2);
    
    GBdist = h5read('/Users/b119user/Downloads/Ajey/Al_Feb2018/20180109_ts3Al_07/GBNormals_v1.1/output/ANGfiles_Alroomtemp.h5','/GBdistance');
    GrainID = h5read('/Users/b119user/Downloads/Ajey/Al_Feb2018/20180109_ts3Al_07/GBNormals_v1.1/output/ANGfiles_Alroomtemp.h5','/GID');
    GBcoord = [newcoord1(GBdist==0),newcoord2(GBdist==0)];
    GBID = GrainID(GBdist==0);
    GBnormals_GB = GBnormals(GBdist==0,:);
    GBTangents_GB = GBTangents(GBdist==0,:);
    
    quiver(GBcoord(:,1),GBcoord(:,2),GBnormals_GB(:,1),GBnormals_GB(:,2))
    
    % gradu = gradient(u,
    % GBTangents = cross(GBnormals,repmat([0,0,1],[length(GBcoord),1]),2);
    figure,
    quiver(GBcoord(:,1),GBcoord(:,2),GBTangents_GB(:,1),GBTangents_GB(:,2))
    mantlecoordidx = rangesearch([newcoord1(:),newcoord2(:)],[GBcoord(:,1),GBcoord(:,2)],50);
    mantlecoordidx_all = [mantlecoordidx{:}].';
    mantlecoord = [newcoord1(mantlecoordidx_all),newcoord2(mantlecoordidx_all)];
    
    mantletangentidx = knnsearch(GBcoord,mantlecoord,'K',1);
    mantleID = GBID(mantletangentidx);
    MantleTangents = [GBTangents_GB(mantletangentidx,1),GBTangents_GB(mantletangentidx,2)];
    MantleNormals = [GBnormals_GB(mantletangentidx,1),GBnormals_GB(mantletangentidx,2)];
    mantleuidx = knnsearch([DICx(:),DICy(:)],[mantlecoord(:,1),mantlecoord(:,2)],'K',1);
    % gradu = gradient(u,0.006);
    % gradv = gradient(v,0.0314);
    mantleu = u_corrected(mantleuidx).*0.0827;
    mantlev = v_corrected(mantleuidx).*0.0827;
    dispmag = dot([mantleu,mantlev],MantleTangents(:,1:2),2);
    
    % figure('pos',[10,10,1800,600])
    % scatter(mantlecoord(:,1),mantlecoord(:,2),[],dispmag,'filled')
    % colormap(redblue)
    % caxis([-200,200])
    % colorbar
    % xlim([0,max(mantlecoord(:,1))])
    % title('GB sliding magnitude','FontSize',20)
    % saveas(gcf,'slidingmag_redblue.png')
    % [IDx5,dist5] = knnsearch(mantlecoord,mantlecoord,'K',3032);
    % diffdisp = zeros(length(mantlecoord),1);
    % for ii = 1:length(mantlecoord)
    %     buff = find(mantleID(IDx5(ii,:)) ~= mantleID(IDx5(ii,1)));
    %     diffdisp(ii) = dispmag(ii) + dispmag(IDx5(ii,buff(1)));
    % end
    % saveas(gcf,'GBslidingmag.png')
    
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
    GBID_mantle = GBID_all(mantletangentidx);
    GBIDpair_mantle = GBID_pair(mantletangentidx,:);
    
    % GBmantleunique = unique(GBIDpair_mantle,'rows');
    % figure,
    % for ii = 1:length(unique(GBIDpair_mantle,'rows'))
    %     GBuniquemantleidx = union(find(all(bsxfun(@eq,GBIDpair_mantle,(GBmantleunique(ii,:))), 2)),find(all(bsxfun(@eq,GBIDpair_mantle,fliplr(GBmantleunique(ii,:))), 2)));
    %     GBuniquecoords{ii} = mantlecoord(GBuniquemantleidx,:);
    %     GBuniquenormal(ii,:) = MantleTangents(GBuniquemantleidx(1),:);
    %     plot(GBuniquecoords{ii}(:,1),GBuniquecoords{ii}(:,2),'.')
    %     hold on
    %     plot(mean(GBuniquecoords{ii}(:,1)),mean(GBuniquecoords{ii}(:,2)),'ro')
    %     GBmidpoint(ii,:) = [mean(GBuniquecoords{ii}(:,1)),mean(GBuniquecoords{ii}(:,2))];
    % end
    %
    % hold off
    for ii = 1:length(MantleTangents)
        if(dot(MantleTangents(ii,:),[1,0]) < 0)
            MantleTangents(ii,:) = -MantleTangents(ii,:);
        end
    end
    for ii = 1:length(MantleNormals)
        if(dot(MantleNormals(ii,:),[1,0]) < 0)
            MantleNormals(ii,:) = -MantleNormals(ii,:);
        end
    end
    
    dispmag = dot([mantleu,mantlev],MantleTangents(:,1:2),2);
    
    figure,
    scatter(mantlecoord(:,1),mantlecoord(:,2),[],dispmag)
    
    GBunique = unique(GBID_pair,'rows');
    %%
    count = 0;
    for ii = 1:length(GBunique)
        
        if (ii>1 && ismember(fliplr(GBunique(ii,:)),GBunique(1:ii-1,:),'rows'))
            continue
        else
            figure('pos',[10,10,1000,800])
            count = count + 1;
            GBuniqueidx = union(find(all(bsxfun(@eq,GBID_pair,(GBunique(ii,:))), 2)),find(all(bsxfun(@eq,GBID_pair,fliplr(GBunique(ii,:))), 2)));
            GBuniquecoords{ii} = GBcoord(GBuniqueidx,:);
            GBuniquenormal(ii,:) = GBnormals_GB(GBuniqueidx(1),:);
            %     plot(GBuniquecoords{ii}(:,1),GBuniquecoords{ii}(:,2),'.k')
            %     hold on
            sortedGBcoords = sortrows(sortrows(GBuniquecoords{ii},2),1);
            pickpointsidx = ceil(linspace(1,length(sortedGBcoords),8));
            pickedpoints = sortedGBcoords(pickpointsidx(3:end-2),:);
            %     plot(pickedpoints(:,1),pickedpoints(:,2),'ro')
            %     GBmidpoint(ii,:) = [mean(GBuniquecoords{ii}(:,1)),mean(GBuniquecoords{ii}(:,2))];
            subplot(3,1,1),
            
            plot(GBcoord(:,1),GBcoord(:,2),'k.')
            %     [pickedpointx, pickedpointy] = getpts;
            hold on
            for jj = 1:length(pickedpoints)
                pickedpointidx = knnsearch(mantlecoord,[pickedpoints(jj,1), pickedpoints(jj,2)]);
                pickedpointsmantle = mantlecoord(pickedpointidx,:);
                if(pickedpointsmantle(1)~=0)
                    slope = MantleNormals(pickedpointidx,2)/MantleNormals(pickedpointidx,1);
                else
                    slope = 1e6;
                end
                linetangent(count,:,jj) = MantleTangents(pickedpointidx,:);
                c = pickedpointsmantle(2) - slope*pickedpointsmantle(1);
                %         kx = linspace(pickedpointsmantle(1)-20,pickedpointsmantle(1)+20,20);
                %         ky = slope.*kx + c;
                t = linspace(-40,40,20);
                for kk = 1:20
                    kx(kk) = pickedpointsmantle(1) + t(kk)*MantleNormals(pickedpointidx,1);
                    ky(kk) = pickedpointsmantle(2) + t(kk)*MantleNormals(pickedpointidx,2);
                end
                % lineidx = knnsearch([DICx(:),DICy(:)],[kx.',ky.']);
                lineidx = knnsearch(mantlecoord,[kx.',ky.']);
                lineu(:,jj) = dispmag(lineidx);
                %     linev = dispmag_subtracted(lineidx);
                reallinecoord = [mantlecoord(lineidx,1),mantlecoord(lineidx,2)];
                plot(reallinecoord(:,1),reallinecoord(:,2),'LineWidth',2.5)
                
                
            end
            titlestr = sprintf('Selected points along GB %d',count);
            title(titlestr,'FontSize',24)
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            
            
            subplot(3,1,2),
            for jj = 1:size(lineu,2)
                plot(linspace(-40,40,20),-(lineu(:,jj)-max(lineu(:,jj))),'o-','LineWidth',5)
                hold on
            end
            plot(linspace(-40,40,20),-(mean(lineu,2)-max(mean(lineu,2))),'o-','LineWidth',5)
            %     line([0,0],[0,max(-(lineu(:,jj)-max(lineu(:,jj))))],'color','k')
            legend('line 1','line 2','line 3','line 4','Mean','location','best')
            xlabel('Distance from GB')
            ylabel('Resolved displacement (\mum)')
            
            xtickvals = linspace(-40,40,5);
            set(gca,'xtick',xtickvals)
            set(gca,'xticklabel',{'-40','-20','GB','20','40'})
            set(gca,'FontSize',20)
            
            subplot(3,1,3),
            scatter(mantlecoord(:,1),mantlecoord(:,2),[],dispmag,'filled')
            hold on
            plot(GBcoord(:,1),GBcoord(:,2),'k.')
            caxis([min(lineu(:)),max(lineu(:))])
            xlim([min(GBuniquecoords{ii}(:,1))-50,max(GBuniquecoords{ii}(:,1))+50])
            ylim([min(GBuniquecoords{ii}(:,2))-50,max(GBuniquecoords{ii}(:,2))+50])
            xticks([])
            yticks([])
            colormap jet
            colorbar            
            
            savefilename = sprintf('GBsliding_line_onlydisp/GBpoints_lines%d_%d.png',count,outerloop);
            saveas(gcf,savefilename)
            close all
        end
    end
    hold off
    
    for jj = 1:length(linetangent)
        GBori(jj,:) = acosd(dot([mean(linetangent(jj,1,:)),mean(linetangent(jj,2,:))],[1,0]));
    end
end
% figure,
% scatter(mantlecoord(:,1),mantlecoord(:,2),[],dispmag)
% colorbar
%%
% figure('pos',[10,10,1600,800]),
% clear lineu
% plot(GBcoord(:,1),GBcoord(:,2),'.')
% [pickedpointx, pickedpointy] = getpts;
% hold on
% for ii = 1:length(pickedpointx)
%     pickedpointidx = knnsearch(mantlecoord,[pickedpointx(ii), pickedpointy(ii)]);
%     pickedpoints = mantlecoord(pickedpointidx,:);
%     if(pickedpoints(1)~=0)
%         slope = MantleNormals(pickedpointidx,2)/MantleNormals(pickedpointidx,1);
%     else
%         slope = 1e6;
%     end
%     c = pickedpoints(2) - slope*pickedpoints(1);
%     kx = linspace(pickedpoints(1)-20,pickedpoints(1)+20,20);
%     ky = slope.*kx + c;
%     % lineidx = knnsearch([DICx(:),DICy(:)],[kx.',ky.']);
%     lineidx = knnsearch(mantlecoord,[kx.',ky.']);
%     lineu(:,ii) = dispmag(lineidx);
% %     linev = dispmag_subtracted(lineidx);
%     reallinecoord = [mantlecoord(lineidx,1),mantlecoord(lineidx,2)];
%     plot(reallinecoord(:,1),reallinecoord(:,2),'.')
%
%
% end
% figure('pos',[10,10,1200,1000]),
% for ii = 1:size(lineu,2)
%     plot(linspace(-20,20,20),-(lineu(:,ii)-max(lineu(:,ii))),'o-')
%     hold on
% end
% legend('line 1','line 2','line 3','line 4','location','best')
% set(gca,'FontSize',20)

toc

