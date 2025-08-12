close all;
I=imread('lena512.bmp');
figure, imshow(I);

I = double(I);
tv = sqrt(diffh(I).^2+diffv(I).^2);
tv_hist
figure, hist(tv(:), 100)

tv_lap = abs(diffh(I))+ abs(diffv(I));
figure, hist(tv_lap(:), 100)

tilefigs;

%%
%pd = fitdist(tv,'Normal')