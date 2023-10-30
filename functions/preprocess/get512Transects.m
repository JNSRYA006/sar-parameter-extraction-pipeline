function [transectData,vertices] = get512Transects(data, topLeft_x, topLeft_y,th, n)
% Inputs:
% data = larger iamge
% topLeft_x = topLeft x coordinate
% topLeft_y = topLeft y coordinate
% th = angle of propogation from horizontal
% n = number of square transects to take
% Outputs:
% transectData =  512 x 512 x n array
% vertices = n x 4 array of starting and ending x and y coordinates

if topLeft_x == 0
    topLeft_x = 1;
    disp("0 input top left x value changed to 1");
end
if topLeft_y == 0
    topLeft_y = 1;
    disp("0 input top left y value changed to 1");
end

transectData = zeros(512,512,n);
vertices = zeros(n,4);

xStart = topLeft_x;
xEnd = topLeft_x + 511;
yStart = topLeft_y;
yEnd = topLeft_y + 511;
angle_incr = 512*(n-1)*tand(th);

for i = 1:n
    transectData(:,:,i) = data(xStart:xEnd,yStart:yEnd);
    vertices(i,1) = xStart;
    vertices(i,2) = yStart;
    vertices(i,3) = xEnd;
    vertices(i,4) = yEnd;
    xStart = xEnd;
    xEnd = xEnd + 511;
    yStart = floor(yStart + angle_incr/(n-1));
    yEnd = yStart + 511;
end
end