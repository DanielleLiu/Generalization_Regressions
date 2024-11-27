function Rsquared = my_Rsquared_coeff(data,data_fit, relativeToMean)
    % R2 correlation coefficient computation
    
    if relativeToMean
        % The total sum of squares
        sum_of_squares = nansum((data-nanmean(data)).^2);
    else %relative to 0, sum of square would be uncentered
            % The total sum of squares
        sum_of_squares = nansum(data.^2);
    end
    
    % The sum of squares of residuals, also called the residual sum of squares:
    sum_of_squares_of_residuals = nansum((data-data_fit).^2);
    
    % definition of the coefficient of correlation is
    Rsquared = 1 - sum_of_squares_of_residuals./sum_of_squares;
    
    
end