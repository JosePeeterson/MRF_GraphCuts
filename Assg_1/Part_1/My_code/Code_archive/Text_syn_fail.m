clear;

% This code is for colour texture from 1 - 6
% texture 1 to 6 is 64x64 & 65x65 pixels so use nW and nH >= 66


nH = 130;
nW = 130;



kS =9;
w = (kS-1)/2;


I = imread('texture11.jpg');



I=double(I)/255;

% missing seedsize, it is the small texture

sigma = kS/6;

G = fspecial('gaussian', [kS kS], sigma); %gaussian weighting kernel
G=G(:); %convert it to 1-d
eps = 0.1; %error threshold

[cH, cW, c] = size(I);

cPatchs=zeros(kS*kS, (cW-kS+1) * (cH-kS+1), c);

for cC = 1:c %for each color channel do
    cPatchs(:,:,cC) = im2col(I(:,:,cC), [kS kS], 'sliding'); %get all possible sliding windows
end

%cPatchs((((kS*kS)+1):3*kS*kS), (cW-kS+1) * (cH-kS+1),:) = 0;

result=zeros(nH,nW,c);

map=false(nH,nW);

result(1:size(I,1),1:size(I,2),:)=I;
    
figure(1)
imshow(result);   

map(1:size(I,1),1:size(I,2),:)=1;

cPatchs=[cPatchs(:,:,1);cPatchs(:,:,2);cPatchs(:,:,3)];

padded_result = padarray(result, [floor(kS/2) floor(kS/2)]);
padded_map = padarray(map, [floor(kS/2) floor(kS/2)]);


G=repmat(G, 1,size(cPatchs, 2)); %repeat it to match the number of candidate patchs
%take into account color channels.
st='G=[';
for cC=1:c-1
    st=strcat(st,'G;');
end
st=strcat(st,'G];');
eval(st);

nex_unfl_row = 1;
nex_unfl_col = 1;

while ~all(all(map))

dilated_map = map;
dilated_map((1:size(I,1)+nex_unfl_row),(1:size(I,2)+nex_unfl_col)) = 1;
surrounding_map = dilated_map - (map);

unfilledIndicies = find(surrounding_map==1);

for p=1:length(unfilledIndicies)


    [py, px] = ind2sub(size(map),unfilledIndicies(p));

    half=floor(kS/2);
    kernel=padded_result(py+floor(kS/2)-half:py+floor(kS/2)+half,px+floor(kS/2)-half:px+floor(kS/2)+half,:);
    mask=padded_map(py+floor(kS/2)-half:py+floor(kS/2)+half,px+floor(kS/2)-half:px+floor(kS/2)+half,:);

    kernel=kernel(:); %1-d
    mask=mask(:);

    st='mask=[';
    for cC=1:c-1
        st=strcat(st,'mask;');
    end
    st=strcat(st,'mask];');
    eval(st);

    kernel=repmat(kernel, 1,size(cPatchs, 2));
    mask=repmat(mask, 1,size(cPatchs, 2));

    dSSD=sum(mask.*G.*(cPatchs-kernel).^2);
    minD = min(dSSD);
    T = minD*(1+eps);

    indicies=find(dSSD<=T);

    index=indicies(randi(length(indicies))); %index
    patch=cPatchs(:,index); %get the patch
    patch=reshape(patch,[kS,kS,c]);

    newPixel = patch(ceil(kS/2),ceil(kS/2),:);
    
    result(py,px,:)=newPixel; %add it to the synthetic texture image
    padded_result(floor(kS/2)+py,floor(kS/2)+px,:)=newPixel; %update the padded result
    map(py,px,:)=1; %raise flag of corresponding pixels in the map
    padded_map(floor(kS/2)+py,floor(kS/2)+px,:)=1;
        
end

nex_unfl_row = nex_unfl_row + 1;
nex_unfl_col = nex_unfl_col + 1;

figure(2)
imshow(result);

end























