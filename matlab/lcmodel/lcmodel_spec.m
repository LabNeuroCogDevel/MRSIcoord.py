function lcmodel_spec(fName)
% LCMODEL_SPECDIR - run lcmodel on fname spectrum output from coordinate placement
% create fname.dir, parse fname, run lcmodel
%
% assume basis function file and lcmodel binary are in the same directory as this script
%
    
%% VARIABLES
lcmodelfolder = fileparts(which(mfilename));
basis = [lcmodelfolder '/' 'basis-sets/gamma_7TJref_te34_297mhz_no2HG_1.basis'];  % for 7 T Jref      
binary_name = [lcmodelfolder '/lcmodel'];
if ~exist(basis, 'file') 
   error(['cannot find basis file: ', basis])
end
if ~exist(binary_name,'file') 
   error(['cannot find lcmodel: ', binary_name])
end

CLEANUP=0;

AX = 297211197;   % for 7 T
% AX = 123151380;     % for 3 T Trio2
% AX = 123254000;     % for 3 T Prisma
ax=num2str(AX);
TE =   17;
TR = 1500; % for Trio
% TR = 2000; % for Prisma
VOLUME = 1;
% NzR = Nslices;
RBW = 3000;
spDT = 333.3000; % (RBW = 1/spDT = 1301 Hz for 3T or 3000 Hz for 7T)
% RBW = 1301; % Trio2
% RBW = 1250; % Prisma
% spDT = 800.0000;    % for 3 T

ppmst = 4.0;
ppmend = 1.8;
LPS = 8; % to enable ps
% LPS = 0;            % suppress ps and coraw
LCORAW = 0;
metab10 = 'mI';
metab9 = 'Gln';
metab8 = 'GABA';
metab7 = 'GSH';
metab6 = 'Glu';
metab5 = 'GPC';
metab4 = 'NAAG';
metab3 = 'NAA';
metab2 = 'Cho';
metab1 = 'Cre';

slsProc = (1:1);
NslsProc = length(slsProc);
% Nkeep = 394; %for ReconRSI_v2
Nkeep = 1024; % for 7T
las = 0;
Jref = 1;
Nx = 1;
Ny = Nx;
NxZ = Nx;
NyZ= NxZ;
IROWST =   1;
IROWEN = NyZ;
ICOLST =   1;
ICOLEN = NxZ;
ISLST  =   1;
ISLEN  = NslsProc;
IPAGE2 = 0;
ISLICE = 1;


rawFname='csi.raw';
controlFname='csi.control';
coordFname='csi.coord';
csvFname = 'spreadsheet.csv';


%%%% ACTUALLY RUN

if ~exist(fName,'file') 
   error(['cannot find inputfile: ', fName, '. (cwd: ',pwd,')'])
end

% dont need to rerun
% TODO: add redo option?
outdir=[fName,'.dir'];
final_out =[outdir '/' csvFname];
if exist(final_out, 'file')
    fprintf("# already have %s. skipping\n", final_out);
    return
end

% try to read data before making directories
dataLC = ifft_spectrum(fName); % also depends on many indexing variable names

prevdir = pwd;
mkdir(outdir);
cd(outdir);

% lcmodel inputs
write_raw(rawFname, dataLC);
write_control(controlFname, fName);
% actually run
system([binary_name ' < ',controlFname]);

% Converting ps to pdf file
if LPS~=0
    %strcmd = ['ps2pdf ' psFname ' ' psFname(1:end-3) '.pdf'];
    %system(strcmd);
    %fprintf('not writting pdf -- ghostscript library error w/matlab\n')
    1;
end

% clean up
if CLEANUP
   delete('*.raw');
   delete('*.control');
end
cd(prevdir)


%%% SUPPORTING FUNCTIONS
function write_raw(fName, dataLC)
    fid=fopen(rawFname,'wt');
    fprintf(fid,' $SEQPAR\n');
    if Jref==1
        fprintf(fid, strcat(' ECHOT = ',num2str(2*TE),'.\n')); % 05/08/15
    else
        fprintf(fid, strcat(' ECHOT = ',num2str(TE),'.\n'));
    end;
    fprintf(fid, strcat(' HZPPPM=',ax(1:3),'.',ax(4:end),'\n'));
    if las==1
        fprintf(fid,' SEQ = ''LASER''\n');
    else
        if Jref==1
            fprintf(fid,' SEQ = ''J-REF''\n');
        else
            fprintf(fid,' SEQ = ''FREQ_REF''\n');
        end;
    end;
    %fprintf(fid, strcat(' ECHOT = ',num2str(TE),'.\n'));
    %fprintf(fid, strcat(' HZPPPM=',ax(1:3),'.',ax(4:end),'\n'));
    %fprintf(fid,' SEQ = ''FREQ_REF''\n');
    fprintf(fid,' $END\n');
    fprintf(fid, strcat(' $NMID ID=''',fName ,''', FMTDAT=''(2e15.6)''\n')); % when running Matlab on Linux
    fprintf(fid, strcat(' TRAMP= 1.0, VOLUME= ',num2str(VOLUME),' $END\n'));
    fprintf(fid,'%15.6E%15.6E\n',dataLC); % when running Matlab on Linux and FMTDAT=''(2e15.6)
    fclose(fid);
end

function write_control(controlFname, fName)
    title = fName;
    if LPS~=0
        psFname='csi.ps';
    end
    if LCORAW~=0
        corawFname='csi.coraw'; 
    end


    fid=fopen(controlFname,'wt');
    fprintf(fid,' $LCMODL\n');
%         fprintf(fid,' OWNER=''Dept of Radiology, Univ of Pittsburgh. Author: Claudiu Schirda, PhD''\n');
%             fprintf(fid,' OWNER=''Department of Radiology, University of Pittsburgh''\n');
% %             fprintf(fid,' Author: ''Claudiu Schirda, PhD - Department of Radiology, University of Pittsburgh''\n');

    fprintf(fid,[' TITLE= ''', title,''' \n']); %%%%
    fprintf(fid,[' FILBAS= ''', basis,'''\n']); %%%%

    fprintf(fid, strcat(' FILRAW=''', rawFname,''' \n'));
    fprintf(fid, strcat(' FILCOO=''', coordFname,''' \n'));
    fprintf(fid, strcat(' FILCSV=''', csvFname,''' \n'));
    if LPS~=0
        fprintf(fid, strcat(' FILPS=''', psFname,''' \n'));
    end
    if LCORAW~=0
        fprintf(fid, strcat(' FILCOR=''', corawFname,''' \n')); 
    end
% % % %             fprintf(fid, strcat(' FILCSV=''', csvFname,''' \n'));
    fprintf(fid, strcat(' HZPPPM=',ax(1:3),'.',ax(4:end),'\n'));

    fprintf(fid,[' NUNFIL=', int2str(Nkeep),'\n']); 
    if RBW<1000
        fprintf(fid,[' DELTAT=.00',int2str(100*spDT),'\n']); % BW<1000 
    else
        fprintf(fid,[' DELTAT=.000',int2str(10*spDT),'\n']); % 1000<BW<10000 
    end

    % TODO: these are variables elsewhere
    fprintf(fid,[' NDCOLS = 1\n']);
    fprintf(fid,[' NDROWS = 1\n']);
    fprintf(fid,[' NDSLIC = 1\n']);
    fprintf(fid,[' ICOLST = 1\n']);
    fprintf(fid,[' ICOLEN = 1\n']);
    fprintf(fid,[' IROWST = 1\n']);
    fprintf(fid,[' IROWEN = 1\n']);
    fprintf(fid,[' ISLICE = 1\n']);

    fprintf(fid,[' neach = 10\n']);  % don't forget 

    fprintf(fid,[' nameac(10) = ''', metab10,''' \n']);
    fprintf(fid,[' nameac(9) = ''', metab9,''' \n']);
    fprintf(fid,[' nameac(8) = ''', metab8,''' \n']);
    fprintf(fid,[' nameac(7) = ''', metab7,''' \n']);
    fprintf(fid,[' nameac(6) = ''', metab6,''' \n']);
    fprintf(fid,[' nameac(5) = ''', metab5,''' \n']);
    fprintf(fid,[' nameac(4) = ''', metab4,''' \n']);
    fprintf(fid,[' nameac(3) = ''', metab3,''' \n']);
    fprintf(fid,[' nameac(2) = ''', metab2,''' \n']);
    fprintf(fid,[' nameac(1) = ''', metab1,''' \n']);


%   fprintf(fid,[' nsimul = 0\n']);  % if don't want MM and Lip
%   fprintf(fid,[' nrefpk = 1\n']);  % for a phantom
%   fprintf(fid,[' ppmref(1,2)= 1.904\n']) 

%   fprintf(fid,[' degzer = 0\n']); % to suppress phase correction
%   fprintf(fid,[' degppm = 0\n']);
%   fprintf(fid,[' sddegz = 0\n']);
%   fprintf(fid,[' sddegp = 0\n']);

    if (LPS~=0 && IPAGE2==0)
        % IPAGE2 = 0 to suppress printing the second page
        %        = 1 (default) to output a second page 
        fprintf(fid,[' IPAGE2 = 0\n']);
    end
    fprintf(fid,[' PPMST=',num2str(ppmst),'\n']);
    fprintf(fid,[' PPMEND=',num2str(ppmend),'\n']);
    fprintf(fid,' LCOORD=9\n');
    if LPS ==-1
        fprintf(fid,[' LPS = 0\n']);
    end
    if LCORAW~=0
        fprintf(fid,' LCORAW=10\n');
    end
        fprintf(fid,' LCSV=11\n');
    fprintf(fid,' $END\n');
    fclose(fid);   
end

function dataLC = ifft_spectrum(fName)
    SP_SLS = zeros(Nkeep,Nx,Ny,NslsProc,'single');
    fid = fopen(fName,'rb');
    SL = fread(fid,'single');
    fclose(fid);
    SL = reshape(SL,2*Nkeep,Nx,Ny);
    SP_SLS(:,:,:,1) = SL(1:Nkeep,:,:) +1i.*SL(Nkeep+1:2*Nkeep,:,:);
    clear SL;
    
    SP_SLS=squeeze(SP_SLS);
    TD_SLS = ifft(SP_SLS,[],1);
    TD_SLS(2:2:end,:,:) = (-1).* TD_SLS(2:2:end,:,:);
    
    % likely TD_SLS(:, 1:1, 1:1, 1:1)
    combDataTDred=TD_SLS(:,IROWST:IROWEN,ICOLST:ICOLEN,ISLST:ISLEN);
    % dataLC_0=reshape(permute(combDataTDred,[1 3 2 4]),2*Nkeep*(ICOLEN-ICOLST+1)*(IROWEN-IROWST+1)*(ISLEN-ISLST+1),1); % 02/17/13 -permute to put the columns (NDCOLS=Nx) in the inner loop, the rows (NDROWS=Ny) in center loop and slices (NDSLIC=NzZ) in the outside loop 
    dataLC_0=reshape(combDataTDred,Nkeep*(ICOLEN-ICOLST+1)*(IROWEN-IROWST+1)*(ISLEN-ISLST+1),1); % 06/13/14 by Victor: do not permute if rows and columns are to correspond to those in sid3  
    dataLC=zeros(2*Nkeep*(ICOLEN-ICOLST+1)*(IROWEN-IROWST+1)*(ISLEN-ISLST+1),1);
    dataLC(1:2:end)=real(dataLC_0);
    dataLC(2:2:end)=imag(dataLC_0);
end

end
