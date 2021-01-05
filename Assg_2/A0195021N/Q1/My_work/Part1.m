% noise_cleaning
clc;
clear;
m_lambda = 10;
fgnd_COLOR = [ 0, 0, 255 ];  % blue foreground 
bck_COLOR = [ 245, 210, 110]; % yellow background

% MRF image with pepper noise
image = imread('Noisy_img.jpg');

% findout the dimensions of image r = height c = width
[r c ~] = size(image);

% create a matrix with replace 0s with empty to save memory.
img_A = sparse(r*c,r*c);

% arrange pixels linearly
lin_pixel = sparse(r*c,2);

% pixel number (total r*c) starting from top left of image
pixel = 1;


for row = 1:r
    for col = 1:c   

    % Find pixel distance between neighbours    
    rgb = [image(row,col,1) image(row,col,2) image(row,col,3)];
    lin_pixel(pixel,1) = dist(fgnd_COLOR,rgb);
    lin_pixel(pixel,2) = dist(bck_COLOR,rgb);
 
      
    neigh_x = col+1;
    neigh_y = row+1;
    img_A(pixel,neigh_x) = m_lambda;
    img_A(neigh_x,pixel) = m_lambda;
    img_A(pixel,neigh_y) = m_lambda;
    img_A(neigh_y,pixel) = m_lambda;
    
      pixel = pixel+1; 
    end
end

% cut graph
[flow,segments] = maxflow(img_A,lin_pixel);
% imshow(res);


reshaped_segments = reshape(segments,c,r)';

% recreate the noise free image
new_image = zeros(r,c,3);
for i = 1:r
    for j = 1:c
        if(reshaped_segments(i,j)==0)
            new_image(i,j,:) = bck_COLOR;
        else
            new_image(i,j,:) = fgnd_COLOR;
        end
    end
end

new_img_rescaled = rescale(new_image);


figure
imshow(new_img_rescaled);
% figure
% imshow(img);

function pxl_dist = dist(c1, c2)
pxl_dist = ( abs( c1(1) - c2(1)) + abs( c1(2) - c2(2) )+ abs( c1(3) - c2(3) )) / 3; 
end