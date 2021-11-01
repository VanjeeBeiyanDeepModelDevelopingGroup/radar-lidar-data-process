function [carte_XSNR,x,y] = polar2carte(polar_xsnr, thetaGrid, rangeGrid, amp)
    [m,n]=size(polar_xsnr);
    thetaGrid = thetaGrid*pi/180;
    [t,r]=meshgrid(thetaGrid,rangeGrid);
%     maxTheta = max(thetaGrid);
%     minTheta = min(thetaGrid);
    maxRange = max(rangeGrid);
    x = linspace(-maxRange/amp,maxRange/amp,n);
    y = linspace(0,maxRange,m);
    [X,Y] = meshgrid(x,y);
    T=atan2(X,Y);
    R=sqrt(X.^2+Y.^2); 
    carte_XSNR = interp2(t,r,polar_xsnr,T,R,'linear',0);
    carte_XSNR = interp2(carte_XSNR,5);
end