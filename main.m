clc;
clear all;close all;
warning('off','all');

%% All global variables
global width;           global height;         
global cell_width;      global cell_height;     
global windowSize;      global SADsize;         
global ref_pixel;       global ref_block;       %Position compensation of pixels and blocks
global searchwindow;    global not_found;       
global average_level;                           
global threshold;                               
global upleft;          global downright;       
global imgpadG;         global imgprepadG;      %The grayscale image of the current image and the previous frame image
global  const1;         global  const2;         

%% 
outputVideo = VideoWriter('_far_result.avi'); 
outputVideo.Quality = 100; 
outputVideo.FrameRate = 5; 
open(outputVideo); 
%% read video 
mov=VideoReader("E:\downloade\SMA-Chip-level-Real-world-Low-Light-Video-Denoising-Based-on-Scalable-Motion-Analysis-main\dataset\noisy_far.avi");
n=mov.NumberOfFrames;
directory = 'E:\downloade\SMA-Chip-level-Real-world-Low-Light-Video-Denoising-Based-on-Scalable-Motion-Analysis-main\speedtest';
if ~exist(directory, 'dir')
    mkdir(directory);
end

%% set save path
ori_pic = fullfile(directory, 'ori_far.jpg');        
denos_pic = fullfile(directory, 'result_far.jpg');    

%% initialize
average_level=5;  
searchwindow=48;  
windowSize=8;  scale=4; 

%total_unmatched = 0; 

ref_pixel=ceil(searchwindow/windowSize+1)*windowSize;
not_found=searchwindow+1;

filter=[1 4 7 4 1;4 16 26 16 4;7 26 41 26 7;4 16 26 16 4;1 4 7 4 1]; 
filter=filter/sum(sum(filter));    

%% processing
for outer_loop=2:n 

    fprintf('The current outer loop variable value is  %d  \n',outer_loop);
    windowSize=8;  
    SADsize=windowSize*2; 
    upleft=(SADsize-windowSize)/2; 
    downright=windowSize+upleft-1; 
    ref_block=ref_pixel/windowSize; 

    
    if outer_loop==2 
        imgpre=read(mov,outer_loop-1); 
        [height,width,~]=size(imgpre); % size

%        cutw=2;cuth=1;
%        width=width/cutw;    
%        height=height/cuth;

        % num of block
        cell_height=floor(height/windowSize); 
        cell_width=floor(width/windowSize);
        height=windowSize*cell_height;
        width=windowSize*cell_width;
        imgpre=imgpre(1:height,1:width,:);% crop

    
        imgprepad=padarray(imgpre,[ref_pixel,ref_pixel],'symmetric','both');% Symmetrical padding
        clear imgpre;
        imgprepadG=double(rgb2gray(imgprepad));% 

        store=uint8(zeros(cell_height*scale,cell_width*scale)); 
        output=cell(cell_height*scale,cell_width*scale);      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %disp(['cell_height: ', num2str(cell_height)]);
        %disp(['cell_width: ', num2str(cell_width)]);
        %disp(['ref_pixel: ', num2str(ref_pixel)]);
        %disp(['windowSize: ', num2str(windowSize)]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mv=zeros(cell_height+ref_pixel/windowSize*2,cell_width+ref_pixel/windowSize*2,2);%Initial mv

    else% 
        cell_height=height/windowSize;
        cell_width=width/windowSize; 
        imgpre=outputD;
        imgprepad=padarray(imgpre,[ref_pixel,ref_pixel],'symmetric','both');  %pad
        %%%%%%%%%%%%%%%%!!!!clear imgpre;
        imgprepadG=imgpadG;    

    end

    %read
    img=read(mov,outer_loop);
    img=img(1:height,1:width,:); 
    imgpad=padarray(img,[ref_pixel,ref_pixel],'symmetric','both');
    %%%%%%%%%%%%clear img;
    imgpadG=double(rgb2gray(imgpad));
    %write    
    %name0=[directory num2str(outer_loop)  ori_pic];
    name0 = fullfile(directory, [num2str(outer_loop) '_ori_far.jpg']);
    %% T
    thresh = getT(imgpadG,imgprepadG,SADsize)*1.5
    %thresh=25.7;
    threshold=thresh*SADsize*SADsize;

    % MEF 
    tic% timing
    mv=mv_initialize(mv);% 

    %MRF
    cell_height=cell_height*2;
    cell_width=cell_width*2; 
    windowSize=windowSize/2;
    SADsize=windowSize*2;
    upleft=(SADsize-windowSize)/2;
    downright=windowSize+upleft-1; 
    ref_block=ref_pixel/windowSize; 
    threshold=thresh*SADsize*SADsize;
    mv_2=mv_refine(mv);

   %full searching
     for k0=0:cell_height*cell_width-1
           i0=floor(k0/cell_width)+1;
           j0=mod(k0,cell_width)+1; 
           i=i0+ref_block;
           j=j0+ref_block;
           flag1=(mv_2(i,j,1)==not_found);
           if flag1
               threshold=thresh*SADsize*SADsize;% reset
               curx=uint16((i-1)*windowSize+1);
               cury=uint16((j-1)*windowSize+1);
               mv_2(i,j,:)=search_T(curx,cury);  
           end
     end    

    %MSF
    cell_height=cell_height*2;
    cell_width=cell_width*2;
    windowSize=windowSize/2;
    SADsize=windowSize*2;
    upleft=(SADsize-windowSize)/2;  
    downright=windowSize+upleft-1;
    ref_block=ref_pixel/windowSize; 
    threshold=thresh*SADsize*SADsize;
    mv_3=mv_refine(mv_2);

    %TF
    count_match=zeros(cell_height,cell_width);
    flag2=(mv_3(ref_block+1:ref_block+cell_height,ref_block+1:ref_block+cell_width,1)~=not_found);
    store=bitset(store,mod((outer_loop-1),average_level)+1,flag2); 
    
    for i=1:average_level
         tmp=double(bitget(store,i));
         count_match=count_match+tmp;
    end 

    count_match=imfilter(count_match,filter,'corr','symmetric','same');
    count_match(count_match<1)=1;

    threshold=thresh*SADsize*SADsize;

  for k0=0:cell_height*cell_width-1 
        i0=floor(k0/cell_width)+1;j0=mod(k0,cell_width)+1;
        i=i0+ref_block;j=j0+ref_block;
        curx=(i-1)*windowSize+1;cury=(j-1)*windowSize+1;
         flag3=(mv_3(i,j,1)~=not_found);
         if flag3
         
            cur_block=imgpad(curx:curx+windowSize-1,cury:cury+windowSize-1,:);          
            mv_block=imgprepad(curx+mv_3(i,j,1):curx+mv_3(i,j,1)+windowSize-1,cury+mv_3(i,j,2):cury+mv_3(i,j,2)+windowSize-1,:);  

            output{i0,j0}=double(cur_block)./count_match(i0,j0)+double(mv_block).*(count_match(i0,j0)-1)./count_match(i0,j0);% temporal.average

         else% match failed
                    sum_inframe=zeros(windowSize,windowSize,3);% 
                    [mv_out2,ct]=search_inframe(curx,cury,thresh*windowSize*windowSize); 
                    mv_out=(mv_out2-1);
                    for k=1:ct
                        sum_inframe=sum_inframe+double(imgpad(mv_out(1,k):mv_out(1,k)+windowSize-1,mv_out(2,k):mv_out(2,k)+windowSize-1,:));
                    end
           output{i0,j0}=sum_inframe./ct; 
         end
  end
    %%%%%%%%%%%%!!!!clear mv_3;clear mv_2;
    outputD=uint8(cell2mat(output));
    name = fullfile(directory, [num2str(outer_loop) '_result_far.jpg']);
    %name=[directory num2str(outer_loop) denos_pic];
    %imwrite(outputD,name);% save
    if outer_loop~=2
       % aviobj=addframe(aviobj,outputD);
       %writeVideo(outputVideo, outputD);  
    end
toc 
end
warning('on','all');
%aviobj=close(aviobj);
close(outputVideo); 



