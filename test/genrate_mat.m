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


%% RegenCoord 
vo = 0;
ho = 0;
angle = 0;
res=216;
% posl = zeros(6, 2); % gui has 6 for left and right
posl(1,1) = 130;
posl(1,2) = 99;
posl(2,1) = 121;
posl(2,2) = 94;
posl = res + 2 - posl
%create transformation matrix, 2rows x 3columns
rotm1 = [cos(angle) sin(angle) vo];
rotm2 = [-sin(angle) cos(angle) ho];

rotm = [rotm1;rotm2]
thirdrow = ones(1,size(posl,1));
posltemp = [posl';thirdrow];
lastrow = [0 0 1];
rotmtemp = [rotm;lastrow];
retpos = rotmtemp\posltemp;
save('test/data/matlab/retpos.mat', 'retpos')

%% Not used - extracted from RegenCoord 
% %generate the rotated and translated scout again
% %this whole section to calculate cropping values
% scoutarray = reshape(scout,1,res*res);
% %point-by-point matrix of loci, 3 rows x res*res columns=basm=scoutmatrx
% basm1 = 1:res;
% basm2 = ones(1,res);
% basm3 = basm2;
% basm = [basm1;basm2;basm3];
% for aa=2:res
%     basm2(1:res)=aa*ones(1,res);
%     basmsub = [basm1;basm2;basm3];
%     basm = [basm,basmsub];
% end
% scoutmatrx=basm;
% xx = rotm*scoutmatrx;  %rotated imagearray (2r x 3c) * (3r x 65536 cols)
% xx = round(xx);         %2rows by 65536cols matrix ea column (row locus,col locus)
% scoutarray2(1:res,1:res)=zeros(res,res);
% %assign the intensities accdg to scout and scoutarray
% for aa=1:res*res
%     if ((xx(1,aa)>=1) & (xx(2,aa)>=1))
%         if ((xx(1,aa)<=256) & (xx(2,aa)<=256))
%         scoutarray2(xx(1,aa),xx(2,aa)) = scoutarray(aa);
%         end
%     end
% end
