function circlemask = ConcentricCircleMaker(data,r1)
%This function makes a series of concentric circles centered at a user
%defined point in a 2 or 3D image.

%Toggle whether or not you want an output image of the concentric circle
PRINT = 0;

%Find the center of concentric circles, select center of nucleus, create
%concentric circle mask
if numel(size(data))>2
    implay(data)
    nucleusframe = input('Specify frame containing center of soma: ');

    [xc,yc] = getpts(imtool(data(:,:,nucleusframe)));

    center = [xc,yc,nucleusframe];
else
    [xc,yc] = getpts(imtool(data));
    center = [xc,yc];
end
imtool close all;
%Determine number of circles required given raius size and image size
if size(data,2)>size(data,1)
    numcircle = ceil(size(data,2)/r1*sqrt(2));
else
    numcircle = ceil(size(data,1)/r1*sqrt(2));
end

%Create 2D circle mask for one image plane. Third dimension refers to which
%annulus we're on
X = ones(size(data,1),1)*(1:size(data,2));
Y = (1:size(data,1))'*ones(1,size(data,2));
Z = (X - xc).^2 + (Y - yc).^2;
circle2 = false(size(data,1),size(data,2));
circle1 = false(size(data,1),size(data,2));
circlemask = false(size(data,1),size(data,2));

for i = 1:numcircle
    circle2(Z <= (r1*i)^2) = 1;
    circle1(Z <= (r1*(i-1))^2) = 1;
    circlemask = circlemask + i*(circle2 - circle1);
    i
end