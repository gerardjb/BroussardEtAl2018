function ColorImg = gray2oneColor(Img,Color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project RGB Image Into One Plane of Color Space
%by: Joey Broussard
%at UC Davis on 07/24/14
%
%This function takes an image in gray color space and converts it to a
%monochromatic version where the color is specified by colorspace input 
%"Color".
%Convert character color specification to numeric representation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2
    Color = [1,0,0];%default color red if not selected
elseif ischar(Color)
    switch Color
        case {'y','yellow'}
            Color = [1 1 0];
        case {'m','magenta'}
            Color = [1 0 1];
        case {'c','cyan'}
            Color = [0 1 1];
        case {'r','red'}
            Color = [1 0 0];
        case {'g','green'}
            Color = [0 1 0];
        case {'b','blue'}
            Color = [0 0 1];
        case {'w','white'}
            Color = [1 1 1];
        case {'k','black'}
            Color = [0 0 0];
        otherwise
            disp('Unrecognized color specifier; using default.');
            Color = [1,0,0];%use default of red if input not recognized
    end
end

%Create Color Space Matrix
R = Color(1)*ones(size(Img));
G = Color(2)*ones(size(Img));
B = Color(3)*ones(size(Img));
%Concatenate color elements
Palette = cat(3,R,G,B);

%Project Image into selected color
for i = 1:3
    ColorImg(:,:,i) = Palette(:,:,i).*double(Img);
end

%Convert to uint8
if max(ColorImg(:)) <= 1 %if coming from mat2gray image, project to 8 bit space
    ColorImg = 255*ColorImg;
end
ColorImg = uint8(ColorImg);