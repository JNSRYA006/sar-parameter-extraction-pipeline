function plotOnMapWind(long,lat,U_10,V_10,sarData)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



%%create a map
worldmap([double(min(lat)) double(max(lat))],[double(min(long)) double(max(long))])
% get the default projection method
mstruct=gcm;
% the map projection method is 'eqdconic'
% check these websites for the project method:
% https://www.mathworks.com/help/map/summary-and-guide-to-projections.html
% https://www.mathworks.com/help/map/eqdconic.html
% By using this projection information, you can convert (lat,lon) to (x,y)
% on the figure
myLat = -32;
myLon = 18;
[x, y] = mfwdtran(mstruct,myLat,myLon);
scatter(x,y,'ro')
hold on
%%image
eLat = -33;
eLong = 17.5;
sizeDeg = 1;
%img = imread(sarData);
latlim=[eLat eLat+sizeDeg];
lonlim=[eLong eLong+sizeDeg];
% show the frame of the image, you can skip this step
[xlim, ylim] = mfwdtran(mstruct,latlim,lonlim);
plot([xlim(1) xlim(2) xlim(2) xlim(1) xlim(1)],[ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)],'y-')
% show the image
R = maprefcells(xlim,ylim,size(sarData));
mapshow(sarData,R)
quiver(long,lat,U_10(:,:,1),V_10(:,:,1),'b','LineWidth',1)
% Vector plot of second time frame
quiver(long,lat,U_10(:,:,2),V_10(:,:,2),'r','LineWidth',1)
hold off
end