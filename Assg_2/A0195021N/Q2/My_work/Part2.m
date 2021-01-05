clc;
clear;
close all

% choose k clusters
k_clusters = 12;

% capture image into RGB
img = im2double(imread('12003.jpg'));

% find dimensions of image array
siz = size(img);

% form K different segments
k_segments = k_clusters;

% choose pixel distance option
pixl_dist = 'sqEuclidean';

% pixel clustering into k regions
pxl_mat = ToMatrix(img); % convert mxnx3 array into mnx3 matrix
% carryout K-means clustering 
[indx, centroid] = kmeans(pxl_mat, k_segments, 'distance', pixl_dist,'maxiter',200);

% calculate data cost for K cluster center
Data_cost = zeros([siz(1:2) k_segments],'single');

for cent_idx = 1:k_segments
    % using covariance matrix per cluster
    icv = inv(cov(pxl_mat(indx==cent_idx,:)));    
    dif_mat = pxl_mat - repmat(centroid(cent_idx,:), [size(pxl_mat,1) 1]);
    % data cost is minus log likelihood of the pixel to belong to each
    % cluster according to its RGB value
    Data_cost(:,:,cent_idx) = reshape(sum((dif_mat*icv).*dif_mat./2,2),siz(1:2));
end

% cut the graph

% smoothness term: 
% constant part
Smooth_cost = ones(k_segments) - eye(k_segments);
% Apply filter convolution parameter
[Hc Vc] = SpatialCues(img);

% cut graph
cut_grph = GraphCut('open', Data_cost, 10*Smooth_cost, exp(-Vc*5), exp(-Hc*5));
[cut_grph L] = GraphCut('expand',cut_grph);
cut_grph = GraphCut('close', cut_grph);

display_k_clusters(L,centroid);

%---------------- Aux Functions ----------------%
function v = ToMatrix(im)
% converts mxnx3 image and returns (mn)x3 matrix
sz = size(im);
v = reshape(im, [prod(sz(1:2)) 3]);
end
%-----------------------------------------------%
function display_k_clusters(L,c)

k_means_clusters = zeros(size(L,1),size(L,2),3);

for i = 1:size(L,1)
    for j = 1:size(L,2)
        index_1 = L(i,j)+1;
        k_means_clusters(i,j,:) = c(index_1,:);
    end
end

figure
imshow(k_means_clusters);

end
%-----------------------------------------------%
%create pre-defined filter.
function [hC vC] = SpatialCues(im)
g = fspecial('gauss', [13 13], sqrt(13));
dy = fspecial('sobel');
vf = conv2(g, dy, 'valid');
sz = size(im);

vC = zeros(sz(1:2));
hC = vC;

for b=1:size(im,3)
    vC = max(vC, abs(imfilter(im(:,:,b), vf, 'symmetric')));
    hC = max(hC, abs(imfilter(im(:,:,b), vf', 'symmetric')));
end

end
