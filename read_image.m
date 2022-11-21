function [I, ni,nj] = read_image(imgPath, isGray, isColor, normalize, scale)

I=double(imread(imgPath));

if isGray && length(size(I)) == 3 && ~isColor
    redChannel = I(:, :, 1);
    greenChannel = I(:, :, 2);
    blueChannel = I(:, :, 3);
    I = .299*double(redChannel) + .587*double(greenChannel) +.114*double(blueChannel);
end

if normalize
    I=I-min(I(:));
    I=I/max(I(:));
end

if scale
    I = imresize(I,[198 NaN]);
end

s = size(I);
ni = s(1);
nj = s(2);

end