function [newcoords1,newcoords2] =registerimages(x1,y1,facedata1,x2,y2,facedata2)
    figure('pos',[10,10,1200,600]),
    surf(x1,y1,facedata1,'LineStyle','None');
    view(2)
    caxis([0,0.02])
    title('Figure 1')
    [FixedPointsx, FixedPointsy] = getpts;
    
    figure('pos',[10,10,1200,600]),
    surf(x2,y2,facedata2,'LineStyle','None');
    view(2)
    title('Figure 2')
    [MovingPointsx, MovingPointsy] = getpts;
    
    fixedPoints = [FixedPointsx, FixedPointsy];
    movingPoints = [MovingPointsx,MovingPointsy];
    tform = fitgeotrans(movingPoints,fixedPoints,'projective');
    [newcoords1,newcoords2] =  transformPointsForward(tform, x2, y2);
    figure,
    plot(x1,y1,'r.')
    hold on
    plot(newcoords1,newcoords2,'g.')
    
    
end