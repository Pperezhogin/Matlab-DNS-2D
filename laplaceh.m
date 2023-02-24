function [ f ] = laplaceh( w )
    f = diff2x(w) + diff2y(w);
end

