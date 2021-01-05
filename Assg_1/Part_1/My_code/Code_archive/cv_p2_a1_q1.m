clear;

texture_image = 'texture1.jpg';

row_multiple = 2;
col_multiple = 2;

small_texture = imread(texture_image);

sm_r = size(small_texture,1);
sm_c = size(small_texture,2);

large_texture = zeros(sm_r*row_multiple,sm_c*col_multiple,3);
large_r = size(large_texture,1);
large_c = size(large_texture,2);

map = zeros(large_r,large_c);

for i = 1:sm_r
    for j = 1:sm_c
        large_texture(i,j,:) = small_texture(i,j,:);
        map(i,j)=1;
    end
end

window_size = ceil(sm_r/4);
gaussian_sigma = window_size/6;

if mod(window_size,2) == 0
    window_size = window_size + 1;
end

maxrow = sm_r;
maxcol = sm_c;

allFilled = 0;
iter = 1;


while(~allFilled)
    
    disp('working')
    
    %choose pixel
    
    max_neighbors = 0;
    chosen_pixel = [1 1];
    chosen_row = 1;
    chosen_col = 1;
    
    for i = 1:(maxrow + 1)
        for j = 1:(maxcol + 1)
            
            if ((i < sm_r && j < sm_c)||(map(i,j)==1))
                continue
            end
            
            % find number of known neighbors
            
            curr_known_neighbors = 0;
            
            for ii = (i - (window_size-1)/2):(i + (window_size - 1)/2)
                for jj = ( j - (window_size - 1)/2) : (j + (window_size - 1)/2)
                    
                    if (ii >= 1 && jj >= 1)
                        if (map(ii,jj)==1)
                            curr_known_neighbors = curr_known_neighbors + 1;
                        end
                    end
                end
            end
            
            if (curr_known_neighbors > max_neighbors)
                max_neighbors = curr_known_neighbors;
                chosen = [i j];
                chosen_row = i;
                chosen_col = j;
                
                if(chosen_row > maxrow)
                    maxrow = chosen_row;
                end
                
                if(chosen_col > maxcol)
                    maxcol = chosen_col;
                end
                
            end
        end
    end
              
                             
    
    % update pixel
    
        % find chosen pixel's neighborhood in the large image, to compare w
        % small image later
        
       chosen_window = zeros(window_size);        
       chosen_window = large_texture((chosen_row - (window_size - 1)/2): (chosen_row + (window_size - 1)/2),(chosen_col - (window_size - 1)/2): (chosen_col + (window_size - 1)/2),:);
       distance_matrix = [];
       gaussian_2D_filter = fspecial('gaussian',window_size,gaussian_sigma);
       
       for i = (1 + (window_size-1)/2):(sm_r - (window_size - 1)/2)
           for j = (1 + (window_size-1)/2):(sm_c - (window_size - 1)/2)
               
               current_st_window = small_texture((i-(window_size - 1)/2):(i+(window_size - 1)/2), (j - (window_size - 1)/2):(j + (window_size - 1)/2), :);
               current_st_window = double(current_st_window);
               current_difference_rgb = gaussian_2D_filter.*(chosen_window - current_st_window);
               
               % for each pixel in the window, find euclidean distance
               % between chosen window and current window
               % store this distance in the distance matrix
               % the distance matrix has 1 column for each current window
               % each column has 1 entry for each pixel in the window
               
%                current_distance_col_vector = zeros(window_size*window_size,1);
               for ii = 1:window_size
                   for jj = 1:window_size
                       
                       dist_ij_vector = [current_difference_rgb(ii,jj,1) current_difference_rgb(ii,jj,2) current_difference_rgb(ii,jj,3)];
                       curr_euclidean_distance_ij = norm(dist_ij_vector);
                       current_distance_col_vector((i-1)*window_size + j) = curr_euclidean_distance_ij;
                       
                       
                   end
               end
               
               distance_matrix = [distance_matrix current_distance_col_vector];
                                             
           end
       end
       
       
       square_distance_matrix = distance_matrix.^2;
       SSD = sum(square_distance_matrix,1);
       
       [minvalue,minindex] = min(SSD);
       
       st_corresp_col = mod(minindex,window_size);
       st_corresp_row = (minindex - st_corresp_col)/window_size+1;
       
%        st_corresp_chosen_window = small_texture((st_corresp_row - (window_size - 1)/2):(st_corresp_row + (window_size - 1)/2),(st_corresp_col - (window_size - 1)/2):(st_corresp_col + (window_size - 1)/2));
       
      

       large_texture(chosen_row,chosen_col,:) = small_texture(st_corresp_row,st_corresp_col,:);
                                        
       % update map
       
       map(chosen_row,chosen_col) = 1;
       
       
       allFilled = all(map);
       disp(iter);
       
       iter = iter + 1;
    
end

figure
imshow(large_texture);