function [topLeftPixelx, topLeftPixely] = latToPixel(latGrid, lonGrid, topLeftLat, topLeftLon)

    [latIndexRow,latIndexCol] = find(round(latGrid,5) == topLeftLat);
    [lonIndexRow,lonIndexCol] = find(round(lonGrid,5) == topLeftLon);
    % Find the indices that are common between latitude and longitude
    commonRowIndices = intersect(latIndexRow, lonIndexRow);
    commonColIndices = intersect(latIndexCol, lonIndexCol);
    
    topLeftPixelx = max(commonRowIndices)-256;
    topLeftPixely = max(commonColIndices)-256;

end