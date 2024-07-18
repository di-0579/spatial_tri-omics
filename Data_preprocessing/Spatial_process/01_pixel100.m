cd './image/P21'


I = imread("./P21_SS.jpg"); I = rgb2gray(I);
BW = imbinarize(I);


figure
imshow(BW)


pixel = 100;
pixel_count = 2*pixel-1;
[numRows,numCols] = size(BW);
pixel_w = numCols/pixel_count;
pixel_h = numRows/pixel_count;
str ="";

for i = 1:100
    x = round(2*(i-1)*pixel_h + 1);
    for j = 1:100
        y = round(2*(j-1)*pixel_w + 1);
        pixel = BW(x:round(x+pixel_w-1),y:round(y+pixel_h-1));
        C = sum(pixel,'all');
        if C > 0
            str = str+','+j+'x'+i;   
        end
    end    
end



fid = fopen('P21_SS.txt','wt');
fprintf(fid, str);
fclose(fid);