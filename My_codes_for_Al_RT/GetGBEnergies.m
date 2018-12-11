clear
close all
tic
% load('strainStep_7.mat');
Eulerangles = squeeze(h5read('/Users/b119user/Downloads/Ajey/Al_Ajey/Test_specimen_output.dream3d','/DataContainers/ImageDataContainer/CellData/EulerAngles'));
Eulerangles = Eulerangles(:,:,:,1);
Eulerangles = permute(Eulerangles,[2,3,1]);

GBdist = h5read('/Users/b119user/Downloads/Ajey/Al_Feb2018/20180109_ts3Al_07/GBNormals_v1.1/output/ANGfiles_Alroomtemp.h5','/GBdistance');
GrainID = h5read('/Users/b119user/Downloads/Ajey/Al_Feb2018/20180109_ts3Al_07/GBNormals_v1.1/output/ANGfiles_Alroomtemp.h5','/GID');
% GBdist(GrainID==23) = NaN;

GBnormals = h5read('/Users/b119user/Downloads/Ajey/Al_Feb2018/20180109_ts3Al_07/GBNormals_v1.1/output/ANGfiles_Alroomtemp.h5','/GBNormals');
% GBnormals = permute(GBnormals,[2,3,1]);
GBnormals = reshape(GBnormals,[3,size(GBdist,1)*size(GBdist,2)]).';
GBnormals = [GBnormals(GBdist==0,1),GBnormals(GBdist==0,2),GBnormals(GBdist==0,3)];
% [newcoord1,newcoord2] = registerimages(nodescoor(:,1),nodescoor(:,2),stress11,x2.',y2.',GrainID.');
% save('registeredcoordinates1.mat','newcoord1')
% save('registeredcoordinates2.mat','newcoord2')
% load('registeredcoordinates1.mat')
% load('registeredcoordinates2.mat')
[y,x] = meshgrid(1:size(GBdist,2),1:size(GBdist,1));
GBcoord = [x(GBdist==0),y(GBdist==0)];

% GBcoord(6151,:) = NaN(1,2);
% GBcoord(6152,:) = NaN(1,2);
% GBcoord(6132,:) = NaN(1,2);

% avgeulang_all = h5read('/Users/b119user/Downloads/Ajey/Al_Feb2018/20180109_ts3Al_07/ANGoutput_corrected','/DataContainers/ImageDataContainer/CellFeatureData/AvgEulerAngles');
% avgeulang_all = avgeulang_all(:,2:end).';
Eulerangles = reshape(Eulerangles,[size(Eulerangles,1)*size(Eulerangles,2),3]);
avgeulangGB = Eulerangles(GBdist==0,:);
% AvgEulGB = avgeulang_all(GBdist==0);

GrainID_GB = GrainID(GBdist==0); 

% GBcoord = [xpos(GBdist==0), ypos(GBdist==0)];
% plot(GBcoord(:,1),GBcoord(:,2),'.')

[IDx5,dist5] = knnsearch(GBcoord,GBcoord,'K',50);

% Mantlecoordidx = rangesearch([CombinedData.x(:), CombinedData.y(:)],GBcoord,3);
% allcoords = [CombinedData.x(:),CombinedData.y(:)];
% exx = CombinedData.exx(:);
% eyy = CombinedData.eyy(:);
% Mantlecoord = [allcoords(unique([Mantlecoordidx{:}]),1),allcoords(unique([Mantlecoordidx{:}]),2)];
% %%

Symm=zeros(3,3,24);
Symm(:,:,1)=[1 0 0; 0 1 0; 0 0 1];
Symm(:,:,2)=[0 0 1; 1 0 0; 0 1 0];
Symm(:,:,3)=[0 1 0; 0 0 1; 1 0 0];
Symm(:,:,4)=[0 -1 0; 0 0 1; -1 0 0];
Symm(:,:,5)=[0 -1 0; 0 0 -1; 1 0 0];
Symm(:,:,6)=[0 1 0; 0 0 -1; -1 0 0];
Symm(:,:,7)=[0 0 -1; 1 0 0; 0 -1 0];
Symm(:,:,8)=[0 0 -1; -1 0 0; 0 1 0];
Symm(:,:,9)=[0 0 1; -1 0 0; 0 -1 0];
Symm(:,:,10)=[-1 0 0; 0 1 0; 0 0 -1];
Symm(:,:,11)=[-1 0 0; 0 -1 0; 0 0 1];
Symm(:,:,12)=[1 0 0; 0 -1 0; 0 0 -1];
Symm(:,:,13)=[0 0 -1; 0 -1 0; -1 0 0];
Symm(:,:,14)=[0 0 1; 0 -1 0; 1 0 0];
Symm(:,:,15)=[0 0 1; 0 1 0; -1 0 0];
Symm(:,:,16)=[0 0 -1; 0 1 0; 1 0 0];
Symm(:,:,17)=[-1 0 0; 0 0 -1; 0 -1 0];
Symm(:,:,18)=[1 0 0; 0 0 -1; 0 1 0];
Symm(:,:,19)=[1 0 0; 0 0 1; 0 -1 0];
Symm(:,:,20)=[-1 0 0; 0 0 1; 0 1 0];
Symm(:,:,21)=[0 -1 0; -1 0 0; 0 0 -1];
Symm(:,:,22)=[0 1 0; -1 0 0; 0 0 1];
Symm(:,:,23)=[0 1 0; 1 0 0; 0 0 -1];
Symm(:,:,24)=[0 -1 0; 1 0 0; 0 0 1];


for ii = 1:length(IDx5)
    
    buff = find(GrainID_GB(IDx5(ii,:)) ~= GrainID_GB(IDx5(ii,1)));
    AvgEulGBall(ii,1:3) = avgeulangGB(ii,:);
    AvgEulGBall(ii,4:6) = avgeulangGB(IDx5(ii,buff(1)),:);
    
    
    A = AvgEulGBall(ii,1:3);
    a = bungeRotationSample2Crystal(A);
    
    %     GBNormals_nonzero(ii,:) = GBNormals(:,NZGBnorm1(ii),NZGBnorm2(ii));
    GBNormals_crystal(ii,:) = a*GBnormals(ii,:).';
    
    B = AvgEulGBall(ii,4:6);
    b = bungeRotationSample2Crystal(B);
    delta_g=b*inv(a);
    for s_op=1:24
        g=Symm(:,:,s_op)*delta_g;
        % Calculate the angle (in radians)
        angle_axis(s_op,1)=acos(0.5*(g(1,1)+g(2,2)+g(3,3)-1));
        % Calculate the axis (Miller indices) upon which the rotation
        % matrix in rotated for lack of a better word by an angle
        % known as angle.
        n_denominator=sqrt((g(2,3)-g(3,2))^2+(g(3,1)-g(1,3))^2+(g(1,2)-g(2,1))^2);
        angle_axis(s_op,2)=(g(2,3)-g(3,2))/n_denominator;
        angle_axis(s_op,3)=(g(3,1)-g(1,3))/n_denominator;
        angle_axis(s_op,4)=(g(1,2)-g(2,1))/n_denominator;
    end
    angle_axis=sortrows(angle_axis,1);
    angle=angle_axis(1,1);
    %     angle = angle*180/pi;
    GB_mis(ii) = angle;
    
    n(1,1)=angle_axis(1,2);
    n(2,1)=angle_axis(1,3);
    n(3,1)=angle_axis(1,4);
    n_normalize=sqrt((n(1,1))^2+(n(2,1))^2+(n(3,1))^2);
    n(1,1)=abs(n(1,1))/n_normalize;
    n(2,1)=abs(n(2,1))/n_normalize;
    n(3,1)=abs(n(3,1))/n_normalize;
    n=sort(n,'ascend');
%     axis(ii,:) = n;
    
     GB_Energy(ii) = GB5DOF((a.'),(b.'),'Al',GBNormals_crystal(ii,:).',-GBNormals_crystal(ii,:).');
%     GB_Energy(ii) = GB5DOF_old(normr(a.'*mat*a),normr(b.'*mat2*b),'Al');
    clear buff
end


fid = fopen('GB_Energy.txt','w');
fprintf(fid,'x\t\ty\tGB Energy\n');
fprintf(fid,'%f %f %f\n',[GBcoord(:,1),GBcoord(:,2),GB_Energy.'].');
fclose(fid);
fid = fopen('GB_misorientation.txt','w');
fprintf(fid,'%f\n',GB_mis.');
fclose(fid);

%%
figure('pos',[10 100 1800 400])
scatter(GBcoord(:,1),GBcoord(:,2),50,GB_Energy,'filled')
axis equal
h = colorbar('westoutside');
axis equal
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xcolor','none')
set(gca,'ycolor','none')
set(gcf,'color','w')
set(gca,'FontSize',20)

title(h,'GB energy (J/m^2)','FontSize',34)
set(gca,'FontSize',34)
colormap jet
% saveas(gcf,'GBenergy.png')

figure('pos',[10 100 1800 400])
scatter(GBcoord(:,1),GBcoord(:,2),50,GB_mis*180/pi,'filled')
axis equal
h = colorbar;
axis equal
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xcolor','none')
set(gca,'ycolor','none')
set(gcf,'color','w')
set(gca,'FontSize',20)
% title(h,'Misorientation angle ($$)','FontSize',34)
title(h,sprintf('Misorientation angle (%c)',char(176)),'FontSize',34)
set(gca,'FontSize',34)
colormap jet
% saveas(gcf,'GBmis.png')

% GBcoordnew_all = intersect();
toc