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


%% K space -- IFFTData.m would write kspace.1.1 (currcvr,curslice)
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

%% make spectrum - ReconCoordinates3.m

%%%%%
%% ENTER HERE 20200410
% -- will need to follow to SpatialTransform2D.m
%%%%%

% in function looks like posltemp comes from
% [posltemp,posrtemp]=RegenCoor(scout,scoutf,posl,posr);
posltemp = retpos

posltemp = res + 2 - posltemp;
posltemp = posltemp';
poslp(1:numrecon,1:2)=posltemp(1:numrecon,1:2)
% convert the offsets to fractional pixel shifts
for m=1:numrecon
    poslpp(m,1) = (poslp(m,1)-1-(res/2))*(rows/res);
    poslpp(m,2) = (poslp(m,2)-1-(res/2))*(cols/res);
end
poslpp = -1*poslpp + 0.5;   %half px shift in both r and c dxns

for m=1:numrecon
    % do the reconstruction
    % convert the data matrix to appropriate formt

    % read kspace -- created above
    outputfilename = CreateName(7); 
    fp2 = fopen(outputfilename,'r');
    SI = fread(fp2,[points*2 rows*cols],'float');
    fclose(fp2);
    SI = SpatialTransform2D(SI);
    ptr = (rows*cols + rows)/2.0 +1;
    spectrum(1:points*2) = SI(1:points*2,ptr);
    %totalshift = (poslpp(m,1)+poslpp(m,2))*3.14159;
    (poslpp(m,1)+poslpp(m,2))
    disp('left, and in 256 res coords');
    poslp(m,:)
    disp('for 240 res coords, rather than 256 ');
    (240/256)*poslp(m,:)
    disp('fraction px shift ')
    poslpp(m,:)
    totalshift = 0;

    srd(1:points) = spectrum(1:points);
    sid(1:points) = spectrum(points+1:points*2);
    scd = srd + i*sid;
    scd = scd * exp(-i*totalshift);
    spectrum(1:points) = real(scd);
    spectrum(points+1:points*2) = imag(scd);

    % PlotSpectrum(spectrum);
    filename = cat(2,dirname,typename,ext);
    fp11 = fopen(filename,'w');
    fwrite(fp11,spectrum,'float');
    fclose(fp11);
end


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
