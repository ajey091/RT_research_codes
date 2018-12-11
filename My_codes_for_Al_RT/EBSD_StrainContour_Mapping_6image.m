clear
clc
close all
%used for FGN5
%% Inputs

%%F11 1 cycle
Strain_Map_File = '_11_oly_200.mat'; % .mat file containing Strain_Data data structure
EBSD_Data_File = '1.ang'; % EBSD Euler angle file
Output_File = 'Al_Sept17_DIC_EBSD'; % Output file name (no extension)

%%
GB_tol = 2.3; % grain boundary misorientation angle tolerance [deg]

                  
%                 % The grid is the smalles rectangle which contains the fiducial marks
FiducialMarks = [   0,  0; % (x, y) coordinates
                    0,  3000; % click in this order, where x and y are the micron measurements
                    15000, 3000;
                    15000, 0];

Nx = 5000; % number of points along x
Ny = 1000; % number of points along y

%
% plotting conventions
%
setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','inToPlane');

N_FiducialMarks = size(FiducialMarks,1);

RangeDown=-0.15;                 %Bottom Range for Colorbar on plot
RangeUp=0.15;                     %Upper Range for Colorbar on plot
Plot='exx';         %What you want to plot. Following are accepted:
%x, y, u, v, sigma, exx, eyy, exy, e1, e2, gamma

%% Load strain data

load(Strain_Map_File);

Strain_Data.x = x;
Strain_Data.y = y;
% Strain_Data.e1 = StitchedData.e1;
% Strain_Data.e2 = StitchedData.e2;
Strain_Data.exx = exx;
Strain_Data.exy = exy;
Strain_Data.eyy = eyy;
% Strain_Data.gamma = StitchedData.gamma;
% Strain_Data.sigma = StitchedData.sigma;
%added 10.13.16
%Strain_Data.u = StitchedData.u;
%Strain_Data.v = StitchedData.v;

Strain_M = size(Strain_Data.x,1);
Strain_N = size(Strain_Data.x,2);

%% Load and process EBSD data
ebsd = loadEBSD_ang(EBSD_Data_File, 'convertEuler2SpatialReferenceFrame');

% Flip x
X_min = min(ebsd.x);
X_max =  max(ebsd.x);
ebsd.x = X_min + (X_min-X_max)/(X_max-X_min)*(ebsd.x - X_max);

% Use MTEX to find grain boundaries and grain IDs
[grains, ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',GB_tol*degree);

grains = smooth(grains);

%
% Convert data in lists to m by n matrices and save this raw data in EBSD_Data structure
%
xx = ebsd.x;
yy = ebsd.y;
IDs = ebsd.grainId;
angles = Euler(ebsd.orientations, 'Bunge'); 
avg_angles = Euler(grains.meanRotation, 'Bunge');

EBSD_Data.x = xx;
EBSD_Data.y = yy;
EBSD_Data.ang = angles;
EBSD_Data.GrainID = IDs;
EBSD_Data.AvgAng = avg_angles(IDs,:);

% % Determine dimensions of EBSD data
% tol = 1e-2; % tolerance to detect changes in y-coordinate
% i = 1;
% while abs(yy(i+1)-yy(i)) < tol
%     i = i + 1;
% end
% EBSD_N = i
% EBSD_M = floor((length(yy)/EBSD_N))
% 
% EBSD_Data.x = zeros(EBSD_M,EBSD_N);
% EBSD_Data.y = zeros(EBSD_M,EBSD_N);
% EBSD_Data.ang = zeros(EBSD_M,EBSD_N,3);
% EBSD_Data.AvgAng = zeros(EBSD_M,EBSD_N,3);
% EBSD_Data.GrainID = zeros(EBSD_M,EBSD_N);
% 
% % Loop over data and form matrices
% k = 1;
% for i=1:EBSD_M
%     for j=1:EBSD_N
%         EBSD_Data.x(i,j) = xx(k);
%         EBSD_Data.y(i,j) = yy(k);
%         EBSD_Data.ang(i,j,:) = angles(k,:);
%         EBSD_Data.GrainID(i,j) = IDs(k);
%         EBSD_Data.AvgAng(i,j,:) = avg_angles(IDs(k),:);
%         k = k + 1;
%     end
% end

EBSD_Data.GBx = grains.boundary.x;
EBSD_Data.GBy = grains.boundary.y;

% Plot data and obtain reference points

% Plot grains using EBSD data 
figure(1)
plot(ebsd,ebsd.orientations.angle./degree)
drawnow
hold on
plot(grains.boundary,'linewidth',1.5)
hold off

% 
% Have user locate each fiducial mark in EBSD data
%
% EBSD_refpts = zeros(N_FiducialMarks,2); % (x,y) locations of fiducial marks in EBSD image
% for i = 1:N_FiducialMarks
%     h = msgbox(sprintf('Place point at fiducial mark %d',i));
%     waitfor(h);
%     h = impoint;
%     Loc = wait(h);
%     EBSD_refpts(i,1) = Loc(1);
%     EBSD_refpts(i,2) = Loc(2);
% end
% EBSD_refpts

EBSD_refpts =    [5.4338    0.4072
    5.3445    0.9259
    0.9172    0.5767
    1.1425    0.1993].*1.0e+03 ;
% Plot the strain map
figure(2)
Fig = (Strain_Data.(Plot) - RangeDown)/(RangeUp - RangeDown);
imshow(Fig(end:-1:1,:))
colormap(jet(128));

% 
% Have user locate each fiducial mark in strain data
%
% StrainMap_refpts = zeros(N_FiducialMarks,2); % (x,y) locations of fiducial marks in strain map image
% for i = 1:N_FiducialMarks
%     h = msgbox(sprintf('Place point at fiducial mark %d',i));
%     waitfor(h);
%     h = impoint;
%     Ind = round(wait(h));
%     StrainMap_refpts(i,1) = Strain_Data.x(Strain_M-Ind(2)+1,Ind(1));
%     StrainMap_refpts(i,2) = Strain_Data.y(Strain_M-Ind(2)+1,Ind(1));
% end
% StrainMap_refpts
StrainMap_refpts =[        5325        4325
        6205       10725
       60040        8095
       57305        3155];
       
%% construct mapping grid
xvec = linspace(min(FiducialMarks(:,1)), max(FiducialMarks(:,1)), Nx);
yvec = linspace(min(FiducialMarks(:,2)), max(FiducialMarks(:,2)), Ny);

MappedData.x = zeros(Ny, Nx);
MappedData.y = zeros(Ny, Nx);

for i = 1:Ny
    for j = 1:Nx
        MappedData.x(i,j) = xvec(j);
        MappedData.y(i,j) = yvec(i);
    end
end

%% Construct affine mappings

%
% Map the specified grid onto the EBSD data
%
Fig1_refpts = EBSD_refpts;
Fig2_refpts = FiducialMarks;


V = zeros(N_FiducialMarks-1,2); % vectors in Figure1
U = zeros(N_FiducialMarks-1,2); % vectors in Figure2

for i = 2:N_FiducialMarks
    V(i-1,:) = Fig1_refpts(i,:) - Fig1_refpts(1,:);
    U(i-1,:) = Fig2_refpts(i,:) - Fig2_refpts(1,:);
end

A = V'*U / (U'*U); % transformation matrix
Ainv = inv(A);

fig1_p0 = Fig1_refpts(1,:)';
fig2_p0 = Fig2_refpts(1,:)';

Transform_EBSD = @(X) fig1_p0 + A*(X-fig2_p0); % affine mapping
InvTransform_EBSD = @(X) fig2_p0 + Ainv*(X-fig1_p0);

%
% Map the specified grid onto the Strain data
%
Fig1_refpts = StrainMap_refpts;
Fig2_refpts = FiducialMarks;


V = zeros(N_FiducialMarks-1,2); % vectors in Figure1
U = zeros(N_FiducialMarks-1,2); % vectors in Figure2

for i = 2:N_FiducialMarks
    V(i-1,:) = Fig1_refpts(i,:) - Fig1_refpts(1,:);
    U(i-1,:) = Fig2_refpts(i,:) - Fig2_refpts(1,:);
end

A = V'*U / (U'*U); % transformation matrix
Ainv = inv(A);

fig1_p0 = Fig1_refpts(1,:)';
fig2_p0 = Fig2_refpts(1,:)';

Transform_StrainMap = @(X) fig1_p0 + A*(X-fig2_p0); % affine mapping
InvTransform_StrainMap = @(X) fig2_p0 + Ainv*(X-fig1_p0);


%% Apply transformation

%
% Map EBSD data
%
MappedData.grainIDs = zeros(Ny,Nx);
MappedData.ang = zeros(Ny,Nx);
MappedData.avgAng = zeros(Ny,Nx);

TransformedX = zeros(Ny,Nx);
TransformedY = zeros(Ny,Nx);

% Find locations in EBSD data domain which correspond to points on the mapping grid
for i = 1:Ny
    for j = 1:Nx
        loc = Transform_EBSD([MappedData.x(i,j); MappedData.y(i,j)]);
        TransformedX(i,j) = loc(1);
        TransformedY(i,j) = loc(2);
    end
end

% perform 2d interpolation to determine values at points on mapping grid
% using the corresponding locations in the EBSD data domain
% MappedData.grainIDs = interp2(EBSD_Data.x, EBSD_Data.y, EBSD_Data.GrainID, TransformedX, TransformedY,'nearest');
% MappedData.avgAng(:,:,1) = interp2(EBSD_Data.x, EBSD_Data.y, EBSD_Data.AvgAng(:,:,1), TransformedX, TransformedY,'nearest');
% MappedData.avgAng(:,:,2) = interp2(EBSD_Data.x, EBSD_Data.y, EBSD_Data.AvgAng(:,:,2), TransformedX, TransformedY,'nearest');
% MappedData.avgAng(:,:,3) = interp2(EBSD_Data.x, EBSD_Data.y, EBSD_Data.AvgAng(:,:,3), TransformedX, TransformedY,'nearest');
% MappedData.ang(:,:,1) = interp2(EBSD_Data.x, EBSD_Data.y, EBSD_Data.ang(:,:,1), TransformedX, TransformedY,'linear');
% MappedData.ang(:,:,2) = interp2(EBSD_Data.x, EBSD_Data.y, EBSD_Data.ang(:,:,2), TransformedX, TransformedY,'linear');
% MappedData.ang(:,:,3) = interp2(EBSD_Data.x, EBSD_Data.y, EBSD_Data.ang(:,:,3), TransformedX, TransformedY,'linear');

MappedData.grainIDs = zeros(Ny,Nx);
MappedData.avgAng = zeros(Ny,Nx,3);
MappedData.ang = zeros(Ny,Nx,3);

F = TriScatteredInterp(EBSD_Data.x, EBSD_Data.y, EBSD_Data.GrainID, 'nearest');
MappedData.grainIDs(:,:) = F(TransformedX, TransformedY);
F = TriScatteredInterp(EBSD_Data.x, EBSD_Data.y, EBSD_Data.AvgAng(:,1), 'nearest');
MappedData.avgAng(:,:,1) = F(TransformedX, TransformedY);
F = TriScatteredInterp(EBSD_Data.x, EBSD_Data.y, EBSD_Data.AvgAng(:,2), 'nearest');
MappedData.avgAng(:,:,2) = F(TransformedX, TransformedY);
F = TriScatteredInterp(EBSD_Data.x, EBSD_Data.y, EBSD_Data.AvgAng(:,3), 'nearest');
MappedData.avgAng(:,:,3) = F(TransformedX, TransformedY);
F = TriScatteredInterp(EBSD_Data.x, EBSD_Data.y, EBSD_Data.ang(:,1), 'linear');
MappedData.ang(:,:,1) = F(TransformedX, TransformedY);
F = TriScatteredInterp(EBSD_Data.x, EBSD_Data.y, EBSD_Data.ang(:,2), 'linear');
MappedData.ang(:,:,2) = F(TransformedX, TransformedY);
F = TriScatteredInterp(EBSD_Data.x, EBSD_Data.y, EBSD_Data.ang(:,3), 'linear');
MappedData.ang(:,:,3) = F(TransformedX, TransformedY);

% Map grain boundary data
MappedData.GBx = zeros(size(EBSD_Data.GBx));
MappedData.GBy = zeros(size(EBSD_Data.GBy));
for k = 1:length(EBSD_Data.GBx)
    loc = InvTransform_EBSD([EBSD_Data.GBx(k); EBSD_Data.GBy(k)]);
    MappedData.GBx(k) = loc(1);
    MappedData.GBy(k) = loc(2);
end

%
% Map strain data
%
MappedData.e1 = zeros(Ny,Nx);
MappedData.e2 = zeros(Ny,Nx);
MappedData.exx = zeros(Ny,Nx);
MappedData.exy = zeros(Ny,Nx);
MappedData.eyy = zeros(Ny,Nx);
MappedData.gamma = zeros(Ny,Nx);
MappedData.sigma = zeros(Ny,Nx);
MappedData.u = zeros(Ny,Nx);
MappedData.v = zeros(Ny,Nx);

% Find locations in Strain data domain which correspond to points on the mapping grid
for i = 1:Ny
    for j = 1:Nx
        loc = Transform_StrainMap([MappedData.x(i,j); MappedData.y(i,j)]);
        TransformedX(i,j) = loc(1);
        TransformedY(i,j) = loc(2);
    end
end

%% perform 2d interpolation to determine values at points on mapping grid
% using the corresponding locations in the Strain data domain
% MappedData.e1 = interp2(Strain_Data.x, Strain_Data.y, Strain_Data.e1, TransformedX, TransformedY,'linear');
% MappedData.e2 = interp2(Strain_Data.x, Strain_Data.y, Strain_Data.e2, TransformedX, TransformedY,'linear');
MappedData.exx = interp2(Strain_Data.x, Strain_Data.y, Strain_Data.exx, TransformedX, TransformedY,'linear');
MappedData.exy = interp2(Strain_Data.x, Strain_Data.y, Strain_Data.exy, TransformedX, TransformedY,'linear');
MappedData.eyy = interp2(Strain_Data.x, Strain_Data.y, Strain_Data.eyy, TransformedX, TransformedY,'linear');
% MappedData.gamma = interp2(Strain_Data.x, Strain_Data.y, Strain_Data.gamma, TransformedX, TransformedY,'linear');
% MappedData.sigma = interp2(Strain_Data.x, Strain_Data.y, Strain_Data.sigma, TransformedX, TransformedY,'linear');
%added 10.13.16 JR
% MappedData.u = interp2(Strain_Data.x, Strain_Data.y, Strain_Data.u, TransformedX, TransformedY,'linear');
% MappedData.v = interp2(Strain_Data.x, Strain_Data.y, Strain_Data.v, TransformedX, TransformedY,'linear');

% figure(6)
% subplot(1,2,1)
% imshow(MappedData.ang(end:-1:1,:,:)./(2*pi))
% title('Mapped EBSD Data')
% axis equal
% drawnow
%
% subplot(1,2,2)
% Fig = (MappedData.(Plot) - RangeDown)/(RangeUp - RangeDown);
% imshow(Fig(end:-1:1,:))
% colormap(jet(128));
% title('Mapped Strain Contours')
% axis equal
% drawnow

%% Crop Data

% Show grains on mapped grid
figure()
%imshow(MappedData.grainIDs(end:-1:1,:)/length(grains))
pcolor(MappedData.x, MappedData.y, MappedData.grainIDs)
shading interp
colormap(lines(128));
hold on
plot(MappedData.GBx, MappedData.GBy, 'k.')
hold off
title('Grain IDs on Mapped Grid')
axis equal
drawnow

% Prompt the user to further crop the data by drawing a rectangle
button = questdlg('Would you like to further crop data?', 'Cropping', 'yes', 'no', 'no');
if strcmp(button,'yes')
    h = msgbox('Drag rectangle over desired region to crop');
    waitfor(h);
    h = imrect;
    rectpos = wait(h);
    
    i_lower = Ny - floor(rectpos(2) + rectpos(4) - 1);
    i_upper = Ny - ceil(rectpos(2));
    j_lower = ceil(rectpos(1));
    j_upper = floor(rectpos(1) + rectpos(3) - 1);
else
    i_lower = 1;
    i_upper = Ny;
    j_lower = 1;
    j_upper = Nx;
end

MM = i_upper - i_lower + 1;
NN = j_upper - j_lower + 1;

%
% Crop the data using the specified rectangle
%
CombinedData.GrainID = MappedData.grainIDs(i_lower:i_upper, j_lower:j_upper);
CombinedData.x = MappedData.x(i_lower:i_upper, j_lower:j_upper);
CombinedData.y = MappedData.y(i_lower:i_upper, j_lower:j_upper);

CombinedData.ang = MappedData.ang(i_lower:i_upper, j_lower:j_upper, :);
CombinedData.avgAng = MappedData.avgAng(i_lower:i_upper, j_lower:j_upper, :);


% CombinedData.e1 = MappedData.e1(i_lower:i_upper, j_lower:j_upper);
% CombinedData.e2 = MappedData.e2(i_lower:i_upper, j_lower:j_upper);
CombinedData.exx = MappedData.exx(i_lower:i_upper, j_lower:j_upper);
CombinedData.exy = MappedData.exy(i_lower:i_upper, j_lower:j_upper);
CombinedData.eyy = MappedData.eyy(i_lower:i_upper, j_lower:j_upper);
% CombinedData.gamma = MappedData.gamma(i_lower:i_upper, j_lower:j_upper);
% CombinedData.sigma = MappedData.sigma(i_lower:i_upper, j_lower:j_upper);
%added 10.13.16 JR
CombinedData.u = MappedData.u(i_lower:i_upper, j_lower:j_upper);
CombinedData.v = MappedData.v(i_lower:i_upper, j_lower:j_upper);

%
% Crop grain boundary data
%
x_lower = xvec(j_lower);
x_upper = xvec(j_upper);
y_lower = yvec(i_lower);
y_upper = yvec(i_upper);


locs = [MappedData.GBx, MappedData.GBy];
test = (locs(:,1)>=x_lower).*(locs(:,1)<=x_upper).*(locs(:,2)>=y_lower).*(locs(:,2)<=y_upper);
indices = find(test);

CombinedData.GBx = MappedData.GBx(indices);
CombinedData.GBy = MappedData.GBy(indices);

clear MappedData

%
% Fix grain IDs for cropped data
%
[SortedIDs, Indices] = sort(CombinedData.GrainID(:));

k = 1;
CombinedData.GrainID(Indices(1)) = 1;
for i = 2:MM*NN
    if SortedIDs(i-1) < SortedIDs(i)
        k = k + 1;
    end
    CombinedData.GrainID(Indices(i)) = k;
end

%
% Show final grains and strain map plotted on the mapping grid
%
figure()
%imshow(CombinedData.GrainID(end:-1:1,:)/k);
pcolor(CombinedData.x, CombinedData.y, CombinedData.GrainID)
shading interp
colormap(lines(k));
hold on
plot(CombinedData.GBx, CombinedData.GBy, 'k.')
hold off
title('Final Mapped Grains')
axis equal
drawnow

figure()
%Fig = (CombinedData.(Plot) - RangeDown)/(RangeUp - RangeDown);
%imshow(Fig(end:-1:1,:));
pcolor(CombinedData.x, CombinedData.y, CombinedData.(Plot))
shading interp
caxis([RangeDown, RangeUp])
colormap(jet(128));
hold on
plot(CombinedData.GBx, CombinedData.GBy, 'k.')
hold off
title('Final Mapped Strain Contours')
axis equal
drawnow

%% Write data to a file

button = questdlg('Would you like to write a file?', 'Write File', 'no', '.txt', '.mat', '.mat');

if strcmp(button,'.txt')
    fid = fopen([Output_File,'.txt'], 'w');
    
    fprintf(fid, 'Nx = %d \t Ny = %d \n', NN, MM);
    fprintf(fid,'GrainID ');
    fprintf(fid,'x \t\t ');
    fprintf(fid,'y \t\t ');
    fprintf(fid,'phi_1 \t\t ');
    fprintf(fid,'phi \t\t ');
    fprintf(fid,'phi_2 \t\t ');
    fprintf(fid,'grainAvgPhi_1 \t ');
    fprintf(fid,'grainAvgPhi \t ');
    fprintf(fid,'grainAvgPhi_2 \t');
%     fprintf(fid,'e_1 \t\t ');
%     fprintf(fid,'e_2 \t\t ');
    fprintf(fid,'e_xx \t\t ');
    fprintf(fid,'e_xy \t\t ');
    fprintf(fid,'e_yy \t\t ');
%     fprintf(fid,'gamma \t\t ');
%     fprintf(fid,'sigma \n');
    %added 10.13.16
    fprintf(fid,'u \t\t ');
    fprintf(fid,'v \t\t');
    for i = 1:MM
        for j = 1:NN
            fprintf(fid,'%d \t ',CombinedData.GrainID(i,j));
            fprintf(fid,'%f \t ',CombinedData.x(i,j));
            fprintf(fid,'%f \t',CombinedData.y(i,j));
            fprintf(fid,'%f \t',CombinedData.ang(i,j,1));
            fprintf(fid,'%f \t',CombinedData.ang(i,j,2));
            fprintf(fid,'%f \t',CombinedData.ang(i,j,3));
            fprintf(fid,'%f \t',CombinedData.avgAng(i,j,1));
            fprintf(fid,'%f \t',CombinedData.avgAng(i,j,2));
            fprintf(fid,'%f \t',CombinedData.avgAng(i,j,3));
%             fprintf(fid,'%f \t',CombinedData.e1(i,j));
%             fprintf(fid,'%f \t',CombinedData.e2(i,j));
            fprintf(fid,'%f \t',CombinedData.exx(i,j));
            fprintf(fid,'%f \t',CombinedData.exy(i,j));
            fprintf(fid,'%f \t',CombinedData.eyy(i,j));
%             fprintf(fid,'%f \t',CombinedData.gamma(i,j));
%             fprintf(fid,'%f \t \n',CombinedData.sigma(i,j));
            %added 10.13.16
%             fprintf(fid,'%f \t',CombinedData.u(i,j));
%             fprintf(fid,'%f \t',CombinedData.v(i,j));
            
            
        end
        if mod(i,5) == 0
            fprintf('\nPrinted %d of %d y-slices\n', i, MM);
        end
    end
    
    fclose(fid);
elseif strcmp(button,'.mat')
    save([Output_File,'.mat'], 'CombinedData', '-mat')
end