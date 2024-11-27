function avgNorm = nanvecnorm(data)
    %if data is a matrix, compute column wise norm and ignore nan.
    %compute sqrt(x1^2 + x2^2 +..+xn^2) / n where n is the number of non-nan
    %elements to take into account different matrix could have different
    %number of non-nan elements.
    avgNorm = data.*data; %first get element wise multiplication
    avgNorm = sum(avgNorm,'omitnan'); %sum of square
    avgNorm = sqrt(avgNorm); %sqrt
    nonNanElements = sum(~isnan(data));
    avgNorm = avgNorm./nonNanElements; %now avg over # of non-nan elements
end