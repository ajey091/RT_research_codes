function [newcoords1,newcoords2] =registerimages_DIC(x1,y1,facedata1,x2,y2,facedata2)
% Function maps the grid of the second set of points onto the first
    figure,
    surf(x1.',y1.',facedata1.','LineStyle','None');
    
    view(2)
    xlim([min(x1(:)),max(x1(:))])
    colormap jet
    caxis([0,0.03])
%     xlim([min(x1(:)),max(x1(:))])
    title('Figure 1')
    [FixedPointsx, FixedPointsy] = getpts;
    
    figure,
    scatter(x2,y2,[],facedata2,'filled');
    
    view(2)
    
    
%     scatter(x2,y2,[],facedata2,'filled');
%     xlim([min(x2),max(x2)])
%     xlim([min(x2)+800 max(x2)-800])
    
%     caxis([0,20])
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