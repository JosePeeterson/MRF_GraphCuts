clc;
clear;

% #################################################################
% USER DEFINED PARAMETERS
%
% load 225x300x3 sunflower images
% ensure images have same pixel size. SxSx3.
ImageA = imread('Image_A.jpg');
ImageB = imread('Image_B.jpg');
ImageB = double(ImageB);
%
% define patch size ptc_sz = 9, should be an odd number
% center of patch is pixel.
ptc_sz = 7;
%
% #################################################################



% window w on left, right, top and bottom of pixel = 4 units
w = (ptc_sz - 1)/2;

if mod(ptc_sz,2)==0 
ptc_sz = ptc_sz + 1;
end

% random assignment of patches from Image A to Image B

size_A = [size(ImageA,1),size(ImageA,2)];
size_B = [size(ImageB,1),size(ImageB,2)];

% Iterations = 5
Itr = 5;

% create window provision for pixels on the boundary
ImageA_NaN = nan(size_A(1)+2*w,size_A(2)+2*w,3);
ImageA_NaN(1+w:size_A(1)+w,1+w:size_A(2)+w,:) = ImageA;

% NFF carries random indices of row and column for image B
% randomly assign corresponding row and column number for each pixel in
% image B. NNF contains two matrices one matrix is random row index for a
% particular i (row) and j (column) pixel in image B. Another matrix is 
% random col. index for a particular i (row) and j (column) pixel in image B.
% i = 225, j = 300 ,
% random row index = NNF(i,j,1), 
% rand col. index = NNF(i,j,2);
% We say ImageA_NaN(i,j,1:3) corresponds to ImageB(NNF(i,j,1),NNF(i,j,2,),1:3)
% Indices [1+w,size_A(1)-w] ensures that ImageB indices match ImageA_NaN
% Indices.

% correction ****************
% 1+w,size_A(1)-w should be 1+w,size_A(1)+w
%

NNF = cat(3,randi([1+w,size_A(1)-w],size_B),randi([1+w,size_A(2)-w],size_B));

% Initialize offsets
% Calculate D for the randomly assigned patches

D = inf(size_A(1),size_A(2),3);

for i = 1:size_A(1)
  for j = 1:size_A(2)
    
    % D is a 9x9 patch matrix initially  
    % For each pixel D is the difference between colour values of
    % patch around that pixel from ImageA_NaN and ImageB with
    % randomised indices due to NNF
    % D is the distance in colour values between two patches of a
    % particular pixel in ImageA_NaN
    
    patch_D = ImageA_NaN(w+i-w:w+i+w, w+j-w:w+j+w, 1:3)...
          - ImageB(NNF(i,j,1)-w:NNF(i,j,1)+w, NNF(i,j,2)-w:NNF(i,j,2)+w, 1:3);
    
    % removes nan numbers and converts patch matrix (ideally 81 elements) into a vector
    patch_D_hue = patch_D(~isnan(patch_D(:,:,1)));  
    patch_D_sat = patch_D(~isnan(patch_D(:,:,2)));  
    patch_D_ill = patch_D(~isnan(patch_D(:,:,3)));  
    
    % sum of squared vector element entries divided by no. of elements for
    % a particular pixel in ImageA_NaN
    
    D(i,j,1) = sum(patch_D_hue.^2)/length(patch_D_hue);  
    D(i,j,2) = sum(patch_D_sat.^2)/length(patch_D_sat); 
    D(i,j,3) = sum(patch_D_ill.^2)/length(patch_D_ill); 
        
  end
end

% propagate the offsets from indices with smaller D and use this offsets
% for its neighbours 


for iteration = 1:Itr

odd_iter = mod(iteration,2)==1;
    
% set scan order limits from Top Right or Bottom Left based on odd/even 
% iteration
if odd_iter % odd
   i_lim = 1:size_A(1); 
   j_lim = 1:size_A(2);
else % even
   i_lim = size_A(1):(-1):1; 
   j_lim = size_A(2):(-1):1;
end

% Find average of the hue, sat and ill to convert D(:,:,3) to D(:,:)
D_avg(:,:) = ( D(:,:,1)+D(:,:,2)+D(:,:,3) )/3;



for i_sc = i_lim
  for j_sc = j_lim

  % propagate from top and left for all the pixels    
  if odd_iter %odd
        % D(ii,jj) contains sum of squared differences for each pixel between
        % ImageA_NaN and ImageB 
        
        % center, top, left
        D_prp(1) = D(i_sc,j_sc);
        D_prp(2) = D(max(1,i_sc-1),j_sc);
        D_prp(3) = D(i_sc,max(1,j_sc-1));
        % ind_min is the index of minimum D from center or top or left pixel.
        [~,ind_min] = min(D_prp);
      
        switch ind_min
            
        % compare the current D(i,j) with its neighbours: top and left pixel
        % If the current D(i,j) is the smallest it will be noticed when we
        % go to the next row or column as D(i,j) will become a neighbour in
        % in the next row or column.
        
        % Top pixel has minimum D compared to Center & Left pixels 
        % we find the offset (stored in NFF) at the top pixel and propagate
        % this offset to the neighbouring pixels in NNF for current i,j.
        
        % propagate from top
        case 2
            % As long as indices in NFF are less than 225 for i and 300 for
            % j. 
            % Note: for i offset value (row index) in NFF(:,:,1) the max
            % value 221 inside NFF(:,:,1) is not considered because of 
            % NNF(i_sc-1,j_sc,1)+1+w<=size_B(1). 
            
            if NNF(i_sc-1,j_sc,1)+1+w<=size_B(1) && NNF(i_sc-1,j_sc,2)+w<=size_B(2)
                % offset value from Top pixel is used for current pixel
                NNF(i_sc,j_sc,:) = NNF(i_sc-1,j_sc,:);
                % only i or (row index) offset value is incremented by 1
                % because current is 1 pixel below so we increment i (row
                % index) offset by 1. 
                NNF(i_sc,j_sc,1) = NNF(i_sc,j_sc,1)+1;
                % calculate again the new Patch_D for this pixel
                new_patch_D = ImageA_NaN(w+i_sc-w:w+i_sc+w, w+j_sc-w:w+j_sc+w, 1:3)...
                       - ImageB(NNF(i_sc,j_sc,1)-w:NNF(i_sc,j_sc,1)+w, NNF(i_sc,j_sc,2)-w:NNF(i_sc,j_sc,2)+w, 1:3);
                
                patch_D_hue = new_patch_D(~isnan(new_patch_D(:,:,1)));  
                patch_D_sat = new_patch_D(~isnan(new_patch_D(:,:,2)));  
                patch_D_ill = new_patch_D(~isnan(new_patch_D(:,:,3)));  

                D(i_sc,j_sc,1) = sum(patch_D_hue.^2)/length(patch_D_hue);  
                D(i_sc,j_sc,2) = sum(patch_D_sat.^2)/length(patch_D_sat); 
                D(i_sc,j_sc,3) = sum(patch_D_ill.^2)/length(patch_D_ill);                                      
            end
            
        % propagate from left
        case 3
            if NNF(i_sc,j_sc-1,1)<=size_B(1) && NNF(i_sc,j_sc-1,2)+1+w<=size_B(2)
                NNF(i_sc,j_sc,:) = NNF(i_sc,j_sc-1,:);
                NNF(i_sc,j_sc,2) = NNF(i_sc,j_sc,2)+1;
                new_patch_D = ImageA_NaN(w+i_sc-w:w+i_sc+w, w+j_sc-w:w+j_sc+w, 1:3)...
                       - ImageB(NNF(i_sc,j_sc,1)-w:NNF(i_sc,j_sc,1)+w, NNF(i_sc,j_sc,2)-w:NNF(i_sc,j_sc,2)+w, 1:3);
                               
                patch_D_hue = new_patch_D(~isnan(new_patch_D(:,:,1)));  
                patch_D_sat = new_patch_D(~isnan(new_patch_D(:,:,2)));  
                patch_D_ill = new_patch_D(~isnan(new_patch_D(:,:,3)));  

                D(i_sc,j_sc,1) = sum(patch_D_hue.^2)/length(patch_D_hue);  
                D(i_sc,j_sc,2) = sum(patch_D_sat.^2)/length(patch_D_sat); 
                D(i_sc,j_sc,3) = sum(patch_D_ill.^2)/length(patch_D_ill); 
                
            end
        end
        
    % propagate from bottom right
    else %even

        D_prp(1) = D(i_sc,j_sc);
        D_prp(2) = D(min(i_sc+1,size_A(1)),j_sc);
        D_prp(3) = D(i_sc,min(j_sc+1,size_A(2)));
        % ind_min is the index of minimum D from center or top or left pixel.
        [~,ind_min] = min(D_prp);
        
        % propagate from bottom
        switch ind_min
        case 2
            if ind_min==2 && NNF(i_sc+1,j_sc,1)-1-w>=1 && NNF(i_sc+1,j_sc,2)-w>=1
                NNF(i_sc,j_sc,:) = NNF(i_sc+1,j_sc,:);
                NNF(i_sc,j_sc,1) = NNF(i_sc,j_sc,1)-1;
                new_patch_D = ImageA_NaN(w+i_sc-w:w+i_sc+w, w+j_sc-w:w+j_sc+w, 1:3)...
                       - ImageB(NNF(i_sc,j_sc,1)-w:NNF(i_sc,j_sc,1)+w, NNF(i_sc,j_sc,2)-w:NNF(i_sc,j_sc,2)+w, 1:3);
                               
                patch_D_hue = new_patch_D(~isnan(new_patch_D(:,:,1)));  
                patch_D_sat = new_patch_D(~isnan(new_patch_D(:,:,2)));  
                patch_D_ill = new_patch_D(~isnan(new_patch_D(:,:,3)));  

                D(i_sc,j_sc,1) = sum(patch_D_hue.^2)/length(patch_D_hue);  
                D(i_sc,j_sc,2) = sum(patch_D_sat.^2)/length(patch_D_sat); 
                D(i_sc,j_sc,3) = sum(patch_D_ill.^2)/length(patch_D_ill); 
                
                
            end

            % propagate from right
        case 3
            if ind_min==3 && NNF(i_sc,j_sc+1,1)-w>=1 && NNF(i_sc,j_sc+1,2)-1-w>=1
            % elseif idx==3 && NNF(ii,jj+1,1)-w>=1 && NNF(ii,jj+1,2)-1-w>=1
                NNF(i_sc,j_sc,:) = NNF(i_sc,j_sc+1,:);
                NNF(i_sc,j_sc,2) = NNF(i_sc,j_sc,2)-1;
                new_patch_D = ImageA_NaN(w+i_sc-w:w+i_sc+w, w+j_sc-w:w+j_sc+w, 1:3)...
                       - ImageB(NNF(i_sc,j_sc,1)-w:NNF(i_sc,j_sc,1)+w, NNF(i_sc,j_sc,2)-w:NNF(i_sc,j_sc,2)+w, 1:3);
                                
                patch_D_hue = new_patch_D(~isnan(new_patch_D(:,:,1)));  
                patch_D_sat = new_patch_D(~isnan(new_patch_D(:,:,2)));  
                patch_D_ill = new_patch_D(~isnan(new_patch_D(:,:,3)));  

                D(i_sc,j_sc,1) = sum(patch_D_hue.^2)/length(patch_D_hue);  
                D(i_sc,j_sc,2) = sum(patch_D_sat.^2)/length(patch_D_sat); 
                D(i_sc,j_sc,3) = sum(patch_D_ill.^2)/length(patch_D_ill);
            end
            
            

        end
  end
  end
end

end

% Reconstruction of image A from image B. 

rec_Img = zeros(size(ImageA));

for i = (1+w):ptc_sz:size(ImageA,1)-w
    for j = (1+w):ptc_sz:size(ImageA,2)-w
     rec_Img(i-w:i+w,j-w:j+w, 1:3) = ImageB(NNF(i,j,1)-w:NNF(i,j,1)+w, NNF(i,j,2)-w:NNF(i,j,2)+w, 1:3);
    end
end

rec_Img = uint8(rec_Img);
figure(1),imshow(rec_Img);



