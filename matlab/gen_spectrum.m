function gen_spectrum(siarray1_f, res, posl, outdir)
% GEN_SPECTRUM - make spectrum files from inputs
% siarray1_f='test/data/siarray.1.1';
% res=216;
% posl =  [30 99; 121 94; 113 89]
% outdir = '/tmp'

    sires=24;
    rows=sires;
    cols=sires;
    sipts=1024;
    sisliceno=1; 

    %% binary float read (little indian)
    % read siarray --  like sid3
    offsetptr=0; % change here if you wanted past the first slice
    fp1 = fopen(siarray1_f,'r');
    if fp1 <= 0, error('cannot open %s (from %s)', siarray1_f, pwd); end
    fseek(fp1,offsetptr,'bof');
    SI = fread(fp1,[2*sipts sires*sires],'float');
    fclose(fp1);


    %% K space -- IFFTData.m would write kspace.1.1 (currcvr,curslice)
    fp1 = fopen(siarray1_f,'r');
    SIData=zeros(rows,cols,sipts);
    np=sipts;
    for a=1:rows; for b=1:cols;
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

    %% save/load is lossy. preform to keep consistant 
    kspaceoutname=fullfile(outdir,'kspace.1.1');
    fp2 = fopen(kspaceoutname,'w');
    for a=1:rows
        for b=1:cols
            fwrite(fp2,kspace(a,b,1:np*2),'float');
        end
    end
    fclose(fp2);

    % DEBUG

    %fp2 = fopen(kspaceoutname,'r');
    %kspace_reread = fread(fp2,[points*2 rows*cols],'float');
    %fclose(fp2);
    % kspace_noread = reshape(permute(kspace,[3,2,1]),2048,24*24);
    % imshow(abs(kspace_noread - kspace_reread)); caxis([0, 10^-5]); colormap(parula);
    % this is done later:
    %  SI = kspace_reread;


    %% RegenCoord 
    vo = 0; ho = 0; angle = 0;
    % posl = zeros(6, 2); % gui has 6 for left and right
    posl_rz = res + 2 - posl
    %create transformation matrix, 2rows x 3columns
    rotm1 = [cos(angle) sin(angle) vo];
    rotm2 = [-sin(angle) cos(angle) ho];

    rotm = [rotm1;rotm2]
    thirdrow = ones(1,size(posl_rz,1));
    posltemp = [posl_rz';thirdrow];
    lastrow = [0 0 1];
    rotmtemp = [rotm;lastrow];
    retpos = rotmtemp\posltemp;
    %save('test/data/matlab/retpos_3.mat', 'retpos')

    %% make spectrum - ReconCoordinates3.m

    % in function looks like posltemp comes from
    posltemp = retpos;
    points=np;

    %% TODO: this is the second time res+2-posl happends? is that okay?
    numrecon = size(posltemp,2);
    posltemp = res + 2 - posltemp;
    posltemp = posltemp';
    poslp(1:numrecon,1:2)=posltemp(1:numrecon,1:2);

    % convert the offsets to fractional pixel shifts
    for m=1:numrecon
        poslpp(m,1) = (poslp(m,1)-1-(res/2))*(rows/res);
        poslpp(m,2) = (poslp(m,2)-1-(res/2))*(cols/res);
    end
    poslpp = -1*poslpp + 0.5;   %half px shift in both r and c dxns
    %save('test/data/matlab/pospp.mat', 'poslpp')
    % this is pospp in python, returned from si.pos_shift()


    %% %% continue on with Spatialtransform2d

    SI_siarray = SI; % not needed. hold onto for interactive debuging
    % read kspace -- created above
    % test against
    % kspaceoutname = '/Volumes/Hera/Projects/7TBrainMech/subjs/11743_20190802/slice_PFC/MRSI_roi/raw//kspace.1.1';
    fp2 = fopen(kspaceoutname,'r');
    SI = fread(fp2,[points*2 rows*cols],'float');
    fclose(fp2);

    % DEBUG
    % SI == permute(kspace,[3,2,1]),2048,24*24)
    % imshow(abs(reshape(permute(kspace,[3,2,1]),2048,24*24) - kspSI)); caxis([0, 10^-5]); colormap(parula);


    totspatial = rows*cols;
    totpts = 2* points;

    % get the data 
    % read kspace
    kspSI = SI;

    % convert to complex data
    data1(1:points,1:totspatial)= kspSI(1:points,1:totspatial) + i*kspSI(points+1:points*2,1:totspatial);
    kspData1 = data1; % NOT NEEDED hold onto for interactive 


    %% SI = SpatialTransform2D(SI);
    %% detour into SpatialTransform2D
    toggleon = 1;
    hanningon = 1;
    rotangle = 0;
    flipvert= 0;
    fliphorz= 0;
    flipslices=0;

    % setup the kspace shifting
    shiftvolume = 1;

    %vertshift = 1.49464; %VertShift
    %horzshift = -1.60098;% HorzShift


    %% resume Reconcoordinates3
    for m=1:numrecon


        vertshift = poslpp(m,1);
        horzshift = poslpp(m,2);
        SHIFTMAT=zeros(rows,cols);
        if (shiftvolume==1)
            for mi=1:rows;for n=1:cols
                    angle = (((mi-1)-(rows/2))*horzshift/rows) + (((n-1)-(cols/2))*vertshift/cols);          
                    angle = angle * 2 * pi;
                    SHIFTMAT(mi,n) = exp(i*angle);          
            end;end
        else
            SHIFTMAT(1:rows,1:cols) = 1+0*i;
        end

        posl_ext = sprintf('%d.%d',posl(m,1),posl(m,2));           
        %save(['test/data/matlab/shiftmat.' posl_ext '.mat'], 'SHIFTMAT')


        % setup matrix and do 2D transform spatially spectral point by spectral point
        for a=1:points

            % reorganize the data
            for b=1:cols
                for c=1:rows
                    ptr = c+((b-1)*rows);
                    data2(c,b)=data1(a,ptr);
                end
            end

            % voxel shift the data 
            data2 = data2.*SHIFTMAT;
            data2 = fft2((data2));

            % rotate images
            if (rotangle~=0)
                data2 = imrotate(data2,rotangle,'nearest','crop');
            end

            % put the data back
            for b=1:cols
                for c=1:rows
                    ptr = c+((b-1)*rows);
                    data1(a,ptr) = (data2(c,b));
                end

            end
        end

        % overrides SI was also kspSI
        SI(1:points,1:totspatial) = real(data1(1:points,1:totspatial));
        SI(points+1:points*2,1:totspatial) = imag(data1(1:points,1:totspatial));
        %save(['test/data/matlab/spatialtransform2d_' posl_ext '.mat'], 'SI')

        % reset data1
        data1 = kspSI(1:points,1:totspatial) + i*kspSI(points+1:points*2,1:totspatial);
        clear data2;
        % do the reconstruction
        % convert the data matrix to appropriate formt
        ptr = (rows*cols + rows)/2.0 +1;
        spectrum(1:points*2) = SI(1:points*2,ptr); % this SI is from sptailTransform
        %totalshift = (poslpp(m,1)+poslpp(m,2))*3.14159;
        totalshift = 0;

        srd(1:points) = spectrum(1:points);
        sid(1:points) = spectrum(points+1:points*2);
        scd = srd + i*sid;
        scd = scd * exp(-i*totalshift);
        spectrum(1:points) = real(scd);
        spectrum(points+1:points*2) = imag(scd);

        %save(['test/data/matlab/spectrum_' posl_ext '.mat'], 'spectrum')
        filename = fullfile(outdir,['spectrum.' posl_ext])
        fp11 = fopen(filename,'w');
        fwrite(fp11,spectrum,'float');
        fclose(fp11);
        % PlotSpectrum(spectrum);
    end
end