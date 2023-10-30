function [transectMSEStruct] = transectMSE(transectData)
    % Get the number of matrices (n)
    dataSize = size(transectData, 3);
    
    % Initialize the structure to store MSE results
    transectMSEStruct = struct('Transects', cell(1, dataSize^2), 'mse', cell(1, dataSize^2));
    index = 1;
    % Calculate MSE for each pair of matrices
    for i = 1:dataSize
        for j = 1:dataSize
            
            % Calculate MSE
            RMSE = sum(rmse(transectData(:,:,i),transectData(:,:,j)));
            MSE = RMSE^2;
            
            % Store results in the structure
            transectMSEStruct(index).Transects = sprintf('Transect %d vs. Transect %d', i, j);
            transectMSEStruct(index).mse = MSE;
            index = index + 1;
        end
    end
end