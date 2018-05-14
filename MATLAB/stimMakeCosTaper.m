% stimMakeCosTaper.m
%
% make cosine-tapered circular window
%
% usage:  [ContrastTaper]=stimMakeCosTaper(imgSize,taperSize);
%
% modified from Aaron's mkTaper, to allow varying degrees of taper
%  taperSize = proportion(0 to 1.0) of image radius which is tapered


function [ContrastTaper]=stimMakeCosTaper(imgSize,taperSize)

mid=imgSize/2;
taper = zeros(imgSize,imgSize);
[x,y]=size(taper);
full_circle_radius = (1-taperSize)*imgSize/2;  %0.8*imgSize/2;
tapered_width = taperSize*imgSize/2;

for i = 1:imgSize
    for j = 1:imgSize
        r = sqrt((i - mid).^2 + (j - mid).^2);
        if r < full_circle_radius
            taper(i,j) = 1;
        end
        if r >= full_circle_radius && r < full_circle_radius + tapered_width
            rx = r-full_circle_radius;
            taper(i,j) = (cos(2*pi*rx/(2*tapered_width))+1)/2;
        end
        if r >= full_circle_radius + tapered_width
            taper(i,j) = 0;
        end
    end
end
ContrastTaper=mat2gray(taper);

