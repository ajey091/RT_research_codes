close all
clear
tic
for ii = 0:11
    loadfilename = sprintf('_%.2d_oly_200.mat',ii);
    load(loadfilename);
    avgstrain(ii+1,1) = mean(exx(~isnan(exx)));
end
   
toc