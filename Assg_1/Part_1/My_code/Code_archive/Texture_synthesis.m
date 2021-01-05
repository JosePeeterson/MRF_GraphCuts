%clc;
clear;

% create a large empty image 'Lar_tex' that is factorx original texture 

sml_tex = imread('texture1.jpg');

factor = 1.2;

sml_sz_r = size(sml_tex,1);
sml_sz_c = size(sml_tex,2);
lar_sz_r = floor((factor*sml_sz_r));
lar_sz_c = floor((factor*sml_sz_c));

lar_tex = zeros(lar_sz_r,lar_sz_c,3);
lar_tex = uint8(lar_tex);

% define a patch of many neighbouring pixels for current pixel
% Make the patch of ODD size, w+1+w for row and col. direction,1 is current
% pixel. patch size = 5; w = 2; 
ptc_sz =5;
w = (ptc_sz-1)/2;

% define gaussian weights
sigma = ptc_sz/6; % design parameter
G = fspecial('gaussian', [ptc_sz ptc_sz], sigma); %gaussian weighting kernel
G = G(:); % convert matrix to 1d  vector with column dominance.

% figure(1)
% imshow(sml_tex);

% Merge small texture to Top Left corner of Big image

for i = 1:sml_sz_r
    for j = 1:sml_sz_c
        lar_tex(i,j,:) = sml_tex(i,j,:);
    end
end

  figure(1)
  imshow(lar_tex);



% padding the sml_tex to scan through all the pixels in small texerure
pad_sml_tex = padarray(sml_tex, [floor(ptc_sz/2) floor(ptc_sz/2)]);
pad_lar_tex = padarray(lar_tex, [floor(ptc_sz/2) floor(ptc_sz/2)]);

% copy the large texture with samll texture image at unpadded locations
pad_lar_tex(1+w:sml_sz_r+w,1+w:sml_sz_c+w,:) = lar_tex(1:sml_sz_r,1:sml_sz_c,:);

% add 1 more extra row and column to ensure scan window due to min_r_nap
% and min_c_map are within the pad_lar_Tex because they are increased by
% 1 pixel.
pad_lar_tex(size(pad_lar_tex,1)+1,size(pad_lar_tex,2)+1,:) = 0;


% convert to double for processing 
pad_lar_tex = double(pad_lar_tex);

%imshow(pad_sml_tex);

% initialize all pixels to be unfilled in map the map of unfilled pixels in large_texgure
map = false(lar_sz_r,lar_sz_c);
% set filled pizels to 1
map(1:sml_sz_r,1:sml_sz_c) = 1;
% initialize minimum_row_map, w + 1 is to account for pading in pad_lar_tex
min_r_map = sml_sz_r + w + 1;
% initialize minimum_column_map, w + 1 is to account for pading in pad_lar_tex
min_c_map = sml_sz_c + w + 1;














  
       % initialize the total number of pixels in small image texture
       % each pixel is associalted with a patch of size ptc_sz
       patchs_vec(ptc_sz*ptc_sz,size(sml_tex,1)*size(sml_tex,2),3) = 0;
       
       % likewise convert all patches from pad_sml_tex into
       % 75xTotalNo.Pixels matrix in small image texture and NOT
       % pad_sml_tex so 64 * 64 = 4096 total pixels. NOT 68*68
       for c = 1:3 
       patchs_vec(:,:,c) = im2col(pad_sml_tex(:,:,c), [ptc_sz ptc_sz], 'sliding'); %get all possible sliding windows
       end     
       
       % convert 3d matrix to 2d with R G B as submatrixes
       patchs_vec = [patchs_vec(:,:,1);patchs_vec(:,:,c);patchs_vec(:,:,c)];
             
       % make as many columns of Gaussian weights as TotalNo.Pixels
       % It has the same number of rows a ptc_sz*ptc_Sz
       G=repmat(G, 1,size(patchs_vec, 2));
       % match the R G B
       G = [G;G;G];

       
while (((min_r_map-w) <= lar_sz_r) || ((min_c_map-w) <= lar_sz_c))       


% check column condition before making 
% if ((min_c_map-w) < lar_sz_c)    
% initialize pizels filled at start of a row, min_c_map
% row pixel rounter intialized to 1 because min_c_map is indexed to pad_lar_tex so +w applies
r_pxl_ctr = 1; 

% vector of all the filled columns
fld_cols = 0;


% create a row of filled pixels 
while r_pxl_ctr <= min_c_map
 
% % find largest number of know neighbouring pixels.

% neighbour size vector for keeping the size of pixels along a row, (min_c_map - w)
neigh_sz_r(1:min_c_map-w) = 0;





      for jj = 1:(min_c_map-w) %% min_c_map-w, because jj starts from 1 in the pad_lar_tex
      
      % intialize the neighbourhood size at particular jj to be zero if it 
      % does not enter the while loop    
      neigh_sz_r(jj) = 0;    
      if (~ismember(jj,fld_cols))
           % find number of pixels in patch of ptc_sz != 0
           find_neigh = pad_lar_tex(min_r_map-w:min_r_map+w,jj:jj+w+w,:);
            % convert RGB info into 1 number because we just want to know if
            % non-zero so instead of doing SSD for RGB this is simpler.
           find_neigh(:,:,1) = find_neigh(:,:,1) + find_neigh(:,:,2) + find_neigh(:,:,3); 
           % find indices of non-negative pixel values in the
           % neighbourhood/patch
           find_neigh_ind = find(find_neigh(:,:,1));
           % find the number/size of non-negative numbers
           neigh_sz_r(jj) = length(find_neigh_ind);         
                        
      end 
      
      end
       % count the filled pixels along this min_r_map
       r_pxl_ctr = r_pxl_ctr + 1;
       
       % neigh_col_max_ind is used as index of center pixel of patch in 
       % pad_lar_tex so needs to be added to w.
       % if using with lar_tex just subtract w so neigh_col_max_ind - w
       [~, neigh_col_max_ind] = max(neigh_sz_r);
       neigh_col_max_ind = neigh_col_max_ind + w;
       
       % update vector of filled column indices
       fld_cols = horzcat(fld_cols, neigh_col_max_ind-w); 
       
       % find the closest neighbour to this pixel with index
       % (min_r_map,neigh_col_max_ind) in pad_lar_tex from pad_sml_tex
       
       % convert neigbourhood/patch into a column vector, neigh_max_vec.
       % first R converts to from 1 - 25, then G 26 - 50, then B 51 - 75 so
       % in total there are 75x1 entries.
       neigh_max_vec = reshape(pad_lar_tex(min_r_map-w:min_r_map+w,neigh_col_max_ind-w:neigh_col_max_ind+w,:)...
           ,[],1);
             
       % make as many columns (TotalNo.Pixels) of neigh_max_vec for
       % subtraction with patches_vec
       neigh_max_vec = repmat(neigh_max_vec, 1,size(patchs_vec, 2));
              
       % find the Gaussian weighted sum of squared differences
       % divide by 255 to normalize the pixel size
       GWSSD = sum(G.*((patchs_vec-neigh_max_vec)/255).^2);
       % pixel with minimum differene with the max neighbourhood pixel
       min_dis_pixel = min(GWSSD);
       T = min_dis_pixel*(1+0.1); %to avoid zero match, take more than 1 minimum
        
       indicies = find(GWSSD<=T);
       
       %pick random patch from the accpeted candidate patchs
       index = indicies(randi(length(indicies))); %index
       Match_patch = patchs_vec(:,index); %get the patch
       Match_patch = reshape(Match_patch,[ptc_sz,ptc_sz,3]); %reshape it
       
       %get the new pixel
       newPixel = Match_patch(ceil(ptc_sz/2),ceil(ptc_sz/2),:);
       newPixel = uint8(newPixel);
       
       lar_tex(min_r_map-w,neigh_col_max_ind-w,:) = newPixel;
             
      
end 
 
%  figure(2)        
%  imshow(lar_tex);



%end












% create a column of filled pixels 



%if ((min_r_map-w) < lar_sz_r)
% column pixel counter intialized to 1
c_pxl_ctr = 1; 

% vector of all the filled rows
fld_rows = 0;


% create a column of filled pixels 
% extra pixel filled in during row formation so min_r_map - 1
  while c_pxl_ctr <= min_r_map-1 
 
% % find largest number of know neighbouring pixels.

% neighbour size vector for keeping the size of pixels along a column, (min_r_map - w - 1)
% extra pixel filled in during row formation so min_r_map - w - 1
  neigh_sz_c(1:min_r_map-w-1) = 0;





      for jj = 1:(min_r_map-w-1) %% min_r_map-w-1, because jj starts from 1 in the pad_lar_tex
      
      % intialize the neighbourhood size at particular jj to be zero if it 
      % does not enter the if statement in order to capture the index.    
      neigh_sz_c(jj) = 0;    
      if (~ismember(jj,fld_rows))
           % find number of pixels in patch of ptc_sz != 0
           find_neigh = pad_lar_tex(jj:jj+w+w,min_c_map-w:min_c_map+w,:);
            % convert RGB info into 1 number because we just want to know if
            % non-zero so instead of doing SSD for RGB this is simpler.
           find_neigh(:,:,1) = find_neigh(:,:,1) + find_neigh(:,:,2) + find_neigh(:,:,3); 
           % find indices of nonzero pixel values in the
           % neighbourhood/patch
           find_neigh_ind = find(find_neigh(:,:,1));
           % find the number/size of non-negative numbers
           neigh_sz_c(jj) = length(find_neigh_ind);         
                        
      end 
      
      end
       % count the filled pixels along this col. min_c_map
       c_pxl_ctr = c_pxl_ctr + 1;
       
       % neigh_col_max_ind is used as index of center pixel of patch in 
       % pad_lar_tex so needs to be added to w.
       % if using with lar_tex just subtract w so neigh_col_max_ind - w
       [~, neigh_row_max_ind] = max(neigh_sz_c);
       neigh_row_max_ind = neigh_row_max_ind + w;
       
       % update vector of filled row indices
       fld_rows = horzcat(fld_rows, neigh_row_max_ind-w); 
       
       % find the closest neighbour to this pixel with index
       % (neigh_row_max_ind,min_c_map,) in pad_lar_tex from pad_sml_tex
       
       % convert neigbourhood/patch into a column vector, neigh_max_vec.
       % first R converts to from 1 - 25, then G 26 - 50, then B 51 - 75 so
       % in total there are 75x1 entries.
       neigh_max_vec = reshape(pad_lar_tex(neigh_row_max_ind-w:neigh_row_max_ind+w,min_c_map-w:min_c_map+w,:)...
           ,[],1);
             
       % make as many columns (TotalNo.Pixels) of neigh_max_vec for
       % subtraction with patches_vec
       neigh_max_vec = repmat(neigh_max_vec, 1,size(patchs_vec, 2));
              
       % find the Gaussian weighted sum of squared differences
       % divide by 255 to normalize the pixel size
       GWSSD = sum(G.*((patchs_vec-neigh_max_vec)/255).^2);
       % pixel with minimum differene with the max neighbourhood pixel
       min_dis_pixel = min(GWSSD);
       T = min_dis_pixel*(1+0.1); %to avoid zero match, take more than 1 minimum
        
       indicies = find(GWSSD<=T);
       
       %pick random patch from the sample image
       index = indicies(randi(length(indicies))); %index
       Match_patch = patchs_vec(:,index); %get the patch
       Match_patch = reshape(Match_patch,[ptc_sz,ptc_sz,3]); %reshape it
       
       %get the new pixel
       newPixel = Match_patch(ceil(ptc_sz/2),ceil(ptc_sz/2),:);
       newPixel = uint8(newPixel);
       
       lar_tex(neigh_row_max_ind-w,min_c_map-w,:) = newPixel;
             
      
  end 
 
%   figure(3)        
%   imshow(lar_tex);

min_r_map = min_r_map + 1;
min_c_map = min_c_map + 1;


%end

end
    figure(2)        
    imshow(lar_tex);























  