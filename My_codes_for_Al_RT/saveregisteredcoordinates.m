
for outerloop = 5:11
    clearvars -except outerloop
    close all
    loadfilename = sprintf('_%2.2d_oly_200.mat',outerloop);
    load(loadfilename)
    x = downsample(downsample(x,4).',4).';
    y = downsample(downsample(y,4).',4).';
    exx = downsample(downsample(exx,4).',4).';
    DICx = x.*0.0827;
    DICy = y.*0.0827;
    EBSDx = squeeze(h5read('/Users/b119user/Downloads/Ajey/Al_Ajey/Aug2017/Test_specimen/ANGfiles/Test_specimen_output_new.dream3d','/DataContainers/ImageDataContainer/CellData/X Position'));
    EBSDy = squeeze(h5read('/Users/b119user/Downloads/Ajey/Al_Ajey/Aug2017/Test_specimen/ANGfiles/Test_specimen_output_new.dream3d','/DataContainers/ImageDataContainer/CellData/Y Position'));
    EBSDFID = squeeze(h5read('/Users/b119user/Downloads/Ajey/Al_Ajey/Aug2017/Test_specimen/ANGfiles/Test_specimen_output_new.dream3d','/DataContainers/ImageDataContainer/CellData/FeatureIds'));   
    EBSDFID = flipud(EBSDFID);
    [newcoord1,newcoord2] = registerimages(DICx,DICy,exx,EBSDx,EBSDy,EBSDFID);
    savefilename1 = sprintf('registeredcoordinates1_%d.mat',outerloop);
    savefilename2 = sprintf('registeredcoordinates2_%d.mat',outerloop);
    save(savefilename1,'newcoord1')
    save(savefilename2,'newcoord2')
    
end