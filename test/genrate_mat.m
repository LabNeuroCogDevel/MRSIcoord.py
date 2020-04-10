%
% `save` precision is not float32
%

% calls are distilled versions of 

%% binary float read (little indian)
% scout as mat (despite the .mat, this file is not a standard matlab mat file)
fname='test/data/mprage_middle.mat';
res=216;
fp1 = fopen(fname,'r');
scout = fread(fp1,[res res],'float');
fclose(fp1);
save('test/data/matlab/scout.mat', 'scout')

% read siarray --  like sid3
siname='test/data/siarray.1.1';
sires=24,sipts=1024, sisliceno=1 
offsetptr=0
fp1 = fopen(siname,'r');
fseek(fp1,offsetptr,'bof');
SI = fread(fp1,[2*sipts sires*sires],'float');
fclose(fp1);
save('test/data/matlab/si.mat', 'SI')


%% K space
fp1 = fopen('test/data/siarray.1.1','r');
SIData=zeros(24,24,1024);
np=1024;
for a=1:24; for b=1:24;
        datar(1:np) = fread(fp1,[np],'float');
        datai(1:np) = fread(fp1,[np],'float');
        SIData(a,b,1:np) = datar(1:np) + i*datai(1:np);
end; end
fclose(fp1);
for a=1:np
    temp(:,:) = SIData(:,:,a);
    temp = fftshift(ifft2(temp));
    tempr = real(temp);
    tempi = imag(temp);
    kspace(:,:,a) = tempr(:,:);
    kspace(:,:,np+a) = tempi(:,:);
end
save('test/data/matlab/kspace.mat', 'kspace')