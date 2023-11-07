function y = layerEval(layers,parcol,x)
%LAYEREVAL evaluates a piecewise function using transport coefficients and 
%material layer values.
%
% LAYEREVAL(layers,parcol,x) is a piecewise function utilising transport
% coefficients for the 1D ADR equation contained in layers.
%
% The input layers contains material layer properties and transport
% coefficient values for the corresponding material layers.
%
% The input parcol indicates which column of the layers variable contains
% the transport coefficient values.
%
% The input x is a single point or array of values.

y = zeros(size(x));

for ii = 1:size(layers,1)-1

    y = y + layers(ii,parcol)*heaviside(x - layers(ii,1)) - layers(ii,parcol)*heaviside(x - layers(ii+1,1));

end

end
