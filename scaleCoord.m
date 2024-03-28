function [xScaled] = scaleCoord(x)
% Does a linear transformation of x using transformation parameters tP or
% transforms to the range [0 .. 1].
%
% Input:
% x  ... [nx1] double, unscaled Values
%
% Output:
% x_scaled ... [nx1] double, scaled Values

% do the linear transformation
xScaled = (x - min(x))/ (max(x) - min(x));
end
