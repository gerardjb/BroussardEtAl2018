function [xCDF,yCDF,sortCDF,h] = plotCDF(x,color)
%This function plots a basic CDF of the input vector x. Outputs include the
%x and y values used in generation of the CDF and a handle for the graph

if nargin==1
   color = 'b'; 
end

[xCDF,sortCDF] = sort(x,'Ascend');
yCDF = ((1:numel(x)) - 1)/(numel(x) - 1);

if nargin==1
    h = plot(xCDF,yCDF);
else
    h = plot(xCDF,yCDF,'Color',color);
end