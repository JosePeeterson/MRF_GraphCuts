clc;
clear;

BK_LoadLib();
noisy_img = imread('Noisy_img.jpg');
noisy_img = double(noisy_img);
%imshow(noisy_img)

FG(1,1,1) = 255;
FG(1,1,2) = 0;
FG(1,1,3) = 0;

BG(1,1,1) = 110;
BG(1,1,2) = 210; 
BG(1,1,3) = 245;

FG_lbl = 1;
BG_lbl = 0;

ht = size(noisy_img,1);
wd = size(noisy_img,2);

%distance of pixel to foreground
noise1 = zeros(size(noisy_img,1),size(noisy_img,2));
for i=1:size(noisy_img,1)
    for j=1:size(noisy_img,2)
noise1(i,j) = norm(reshape(noisy_img(i,j,:),[3,1]) - reshape(FG(1,1,:),[3,1]));
    end
end
noise1 = noise1./max(noise1);


%distance of pixel to Background
noise2 = zeros(size(noisy_img,1),size(noisy_img,2));
for i=1:size(noisy_img,1)
    for j=1:size(noisy_img,2)
noise2(i,j) = norm(reshape(noisy_img(i,j,:),[3,1]) - reshape(BG(1,1,:),[3,1]));
    end
end
noise2 = noise2./max(noise2);

noise = [noise1(:)';noise2(:)'];

nb = sparse(wd*ht,wd*ht);

% initialize all nb pixel to pixel edge connections to be foreground points.
for y=1:ht % set up a grid-like neighbourhood, arbitrarily
    for x=1:wd
        if (x < wd), nb((y-1)*wd+x,(y-1)*wd+x+1) = 1; end
        if (y < ht), nb((y-1)*wd+x, y   *wd+x  ) = 1; end
    end
end


distmap_rt = zeros(ht,wd);
distmap_bt = zeros(ht,wd);

for i=1:ht-1
    for j=1:wd-1
    distmap_rt(i,j) = norm(reshape(noisy_img(i,j,:),[3,1])...
                    -reshape(noisy_img(i,j+1,:),[3,1]));
    distmap_bt(i,j) = norm(reshape(noisy_img(i,j,:),[3,1])...
                    -reshape(noisy_img(i+1,j,:),[3,1]));
    end
end
distmap_rt = distmap_rt / max(distmap_rt(:));
distmap_bt = distmap_bt / max(distmap_bt(:));
distmap_rt(i+1,:) = 1;
distmap_bt(:,j+1) = 1;


distmap = [distmap_rt(:)'; distmap_bt(:)'];
%distmap = uint8(distmap);

hinc = BK_Create(wd*ht,(wd*(ht-1) + (wd-1)*ht));
BK_SetNeighbors(hinc,nb);

time_inc = [];
time_new = [];

figure;

lambda = 16;


while (lambda >= 1)

    
    newdc = double(noise+lambda*distmap);
     
    
    BK_SetUnary(hinc,newdc); 
    tic; 
        e_inc = BK_Minimize(hinc);
    time_inc(end+1) = toc;
    lab = BK_GetLabeling(hinc);
    
    imagesc(reshape(lab,[wd,ht]));
    drawnow;
    
    hnew = BK_Create(wd*ht,2*wd*ht);
    BK_SetNeighbors(hnew,nb);
    BK_SetUnary(hnew,newdc); 
    tic; 
        e_new = BK_Minimize(hnew);
    time_new(end+1) = toc;
    BK_Delete(hnew);
    
    Assert(abs(e_inc - e_new) < 1e-6);
    lambda = lambda*0.9;
    
   
end

BK_Delete(hinc);




