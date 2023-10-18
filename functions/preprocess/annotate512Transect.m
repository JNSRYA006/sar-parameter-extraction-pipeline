function [rect,textLabel] = annotate512Transect(x,y,number,colour,text_colour,text_bg)
% annotate512Transect returns a rectangle and text object to annotate transects in a larger image
% The x and y values are the coordinates of the upper left corner of the
% transect
% The number represents the text value to label the transect
% colour is the colour of the square
% text_bg is either 1 or 0 to determine the prescence or absence of a white
% text background
    rect = drawRect(x,y,colour);
    textLabel = addNumberLabel(x, y, number, text_colour, text_bg);
end

function [rect] = drawRect(x,y,colour)
% drawRect returns a rectangle object at the location of the transect
% The x and y values are the coordinates of the upper left corner of the
% transect
    rect = rectangle('Position',[x,y,512,512],'EdgeColor', colour,'LineWidth', 1,'LineStyle','-');
end

function [textLabel] = addNumberLabel(xOfRect, yOfRect, number, colour, bg)
% addNumberLabel returns a text object in the middle of the transect
% The x and y values are the coordinates of the upper left corner of the
% transect
% The number represents the text value to label the transect
% bg is either 1 or 0, where 1 represents the prescence of a white
% background, and 0 represents no background for the text
    if bg
        textLabel = text(xOfRect+256,yOfRect+256,num2str(number),'FontSize',16,'Color',colour,'BackgroundColor','w');
    end
    textLabel = text(xOfRect+256,yOfRect+256,num2str(number),'FontSize',16,'Color',colour);
end