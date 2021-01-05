clc;
clear;

% #############################################################
%  User Defined Parameters
%
% define neighbourchood size or patch size
% for minimum performance please use patch size of 9 and above
ptc_sz = 3;
%
%
% For texture 1 to 6 use scale factor of larger texture 1.5
% For texture 7 to 11 use scale factor of larger texture 1.25 to speed up
% texture formation
scale_factor = 1.5;
%
% define gaussian sigma/width which is standard deviation
sigma = ptc_sz/6;
%
%error threshold to capture at least one error
eps = 0.1; 
%
% capture the small texture image 
sml_tex = imread('texture2.jpg');
%
% #############################################################

% check if patch size is odd, make it odd
if mod(ptc_sz,2)==0
ptc_sz = ptc_sz +1;
end

tic; % start the time measurement for particular neighbourhood size 

% Define window for large texture
lar_tex_H = round(size(sml_tex,1)*scale_factor);
lar_tex_W = round(size(sml_tex,2)*scale_factor);

% Normalize the pixel values to work with small numbers
sml_tex=double(sml_tex)/255;

% TEST
figure(1)
imshow(sml_tex);

G = fspecial('gaussian', [ptc_sz ptc_sz], sigma); %gaussian weighting kernel
G=G(:); %convert it to 1-d

% get the size of small texture, sample patches heigh width and colour
[smpl_H_pat, smpl_W_pat, c] = size(sml_tex);

% sample patch from small texture with pixel in center 
smpl_Patchs=zeros(ptc_sz*ptc_sz, (smpl_W_pat-ptc_sz+1) * (smpl_H_pat-ptc_sz+1), c);

%for each color channel capture the patches into a column vector
for c_Ctr = 1:c 
    smpl_Patchs(:,:,c_Ctr) = im2col(sml_tex(:,:,c_Ctr), [ptc_sz ptc_sz], 'sliding');
end

% display final large texture using this 
% it fills/stores the unfilled pixels
fill_lar_tex = zeros(lar_tex_H,lar_tex_W,c);

% create large texture map of logical 0's for intialization
lar_tex_map=false(lar_tex_H,lar_tex_W);

% get the small texture image into fill_lar_tex
fill_lar_tex(1:size(sml_tex,1),1:size(sml_tex,2),:) = sml_tex;

% TEST
figure(2)
imshow(fill_lar_tex);   

% Inititalize large texture map of logical 1's at filled pixels
lar_tex_map(1:size(sml_tex,1),1:size(sml_tex,2),:) = 1;

% check if image is greyscale or colour to define correct 3d reduction to
% mxn 2d
if c == 3
smpl_Patchs = [smpl_Patchs(:,:,1);smpl_Patchs(:,:,2);smpl_Patchs(:,:,3)];
elseif c == 1
smpl_Patchs = smpl_Patchs(:,:,1);    
end

% pad the fill_lar_tex so that scanning can be done on all the filled pixles
padded_fill_lar_tex = padarray(fill_lar_tex, [floor(ptc_sz/2) floor(ptc_sz/2)]);
padded_lar_tex_map = padarray(lar_tex_map, [floor(ptc_sz/2) floor(ptc_sz/2)]);

%repeat it to match the number of smpl_Patchs
G=repmat(G, 1,size(smpl_Patchs, 2)); 
%take into account color channels.
st='G=[';
for c_Ctr=1:c-1
    st=strcat(st,'G;');
end
st=strcat(st,'G];');
eval(st);

% next unfilled row counter
nex_unfl_row = 1;
% next unfilled column counter
nex_unfl_col = 1;

% iterate through until all the unfilled pixels are filled

while ~all(all(lar_tex_map))

% greatest neighbour/patch map of filled and first to be filled pixels
grt_ptc_map = lar_tex_map;
grt_ptc_map((1:size(sml_tex,1)+nex_unfl_row),(1:size(sml_tex,2)+nex_unfl_col)) = 1;

% map of nearest unfilled pixels with greatest neighbours
nrst_unfld_map = grt_ptc_map - (lar_tex_map);

% get the indices of unfilled pixels in order of greatest neighbours
unfld_Inds = find(nrst_unfld_map==1);

for pxl=1:length(unfld_Inds)

% convert column number to x and y coordinates in ptc_Sz ptc_Sz matrix
    [pxl_y, pxl_x] = ind2sub(size(lar_tex_map),unfld_Inds(pxl));

% width of patch size    
    w=floor(ptc_sz/2);
% window of     
    window = padded_fill_lar_tex(pxl_y+floor(ptc_sz/2)-w:pxl_y+floor(ptc_sz/2)+w,pxl_x+floor(ptc_sz/2)-w:pxl_x+floor(ptc_sz/2)+w,:);
    mask=padded_lar_tex_map(pxl_y+floor(ptc_sz/2)-w:pxl_y+floor(ptc_sz/2)+w,pxl_x+floor(ptc_sz/2)-w:pxl_x+floor(ptc_sz/2)+w,:);

    window = window(:); %1-d
    mask=mask(:); % 1d for GWSSD calculation

    st='mask=[';
    for c_Ctr=1:c-1
        st=strcat(st,'mask;');
    end
    st=strcat(st,'mask];');
    eval(st);

    window=repmat(window, 1,size(smpl_Patchs, 2));
    mask=repmat(mask, 1,size(smpl_Patchs, 2));

    dSSD=sum(mask.*G.*(smpl_Patchs-window).^2);
    minD = min(dSSD);
    T = minD*(1+eps);

    indicies=find(dSSD<=T);

    index=indicies(randi(length(indicies))); %index
    patch=smpl_Patchs(:,index); %get the patch
    patch=reshape(patch,[ptc_sz,ptc_sz,c]);

    newPixel = patch(ceil(ptc_sz/2),ceil(ptc_sz/2),:);
    
    fill_lar_tex(pxl_y,pxl_x,:)=newPixel; %add it to the synthetic texture image
    padded_fill_lar_tex(floor(ptc_sz/2)+pxl_y,floor(ptc_sz/2)+pxl_x,:) = newPixel; %update the padded result
    lar_tex_map(pxl_y,pxl_x,:)=1; %raise flag of corresponding pixels in the map
    padded_lar_tex_map(floor(ptc_sz/2)+pxl_y,floor(ptc_sz/2)+pxl_x,:)=1;
        
end

% update unfilled columns and rows to next row and column 
nex_unfl_row = nex_unfl_row + 1;
nex_unfl_col = nex_unfl_col + 1;

end

%TEST
figure(3)
imshow(fill_lar_tex);


















