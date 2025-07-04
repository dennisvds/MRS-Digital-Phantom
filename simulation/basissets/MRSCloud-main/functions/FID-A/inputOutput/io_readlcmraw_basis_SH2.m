% io_readlcmraw_basis.m
% Jamie Near, McGill University 2014.
%
% USAGE:
% out=io_readlcmraw_basis(filename);
%
% DESCRIPTION:
% Reads entire LCModel .basis file into multiple FID-A data structures in MATLAB.
%
% INPUTS:
% filename   = filename of LCModel .basis file.
%
% OUTPUTS:
% out        = Input basis set saved as a structure in which each field is
%               an individual metabolite basis spectrum in FID-A structure 
%               format.

function [out]=io_readlcmraw_basis_SH(filename)

% Begin to read the data.
fid=fopen(filename);
linenum=1;
line=fgets(fid);

%look for FWHMBA
fwhmba_index=findstr(line,'FWHMBA=');
while isempty(fwhmba_index);
    line=fgets(fid);
    fwhmba_index=findstr(line,'FWHMBA=');
end
linewidth=str2num(line(fwhmba_index+7:end-2));

%look for HZPPM
hzpppm_index=findstr(line,'HZPPPM=');
while isempty(hzpppm_index);
    line=fgets(fid);
    hzpppm_index=findstr(line,'HZPPPM=');
end
hzpppm=str2num(line(hzpppm_index+7:end-2));
Bo=hzpppm/42.577;
linewidth=linewidth*hzpppm;

%look for TE
te_index=findstr(line,'ECHOT=');
while isempty(te_index);
    line=fgets(fid);
    te_index=findstr(line,'ECHOT=');
end
te=str2num(line(te_index+6:end-2));

%look for spectral width
badelt_index=findstr(line,'BADELT=');
while isempty(badelt_index);
    line=fgets(fid);
    badelt_index=findstr(line,'BADELT=');
end
dwelltime=str2num(line(badelt_index+7:end-2));
spectralwidth=1/dwelltime;

fileEnd=false;

while ~feof(fid)
    %Look for the metabolite name
    metab_index=findstr(line,'METABO=');
    while isempty(metab_index);
        line=fgets(fid);
        metab_index=findstr(line,'METABO=');
    end
    %metabName=line(metab_index+10:end-3);
    metabName=line(metab_index+8:end-3);
    
    hdrEnd_index=findstr(line,'$END');
    while isempty(hdrEnd_index);
        line=fgets(fid);
        hdrEnd_index=findstr(line,'$END');
    end
    
    line=fgets(fid);
    
    basisused_index=findstr(line,'$BASIS');
    linenum=1;
    RF=[];

    % If the line is empty skip it
    while ~isempty(line) && isempty(basisused_index) && ~fileEnd
        %dataline=line(1:semicol_index-2);
        [A,count, errmsg, nextindex] = sscanf(line, '%f', inf);
        % If read failed, output the error
        if ~isempty(errmsg);
            fclose(fid);
            line %scnh
            error('READLCMRAW_BASIS failed with read error: %s', errmsg);
        end
        % Store the read values into rf array
        RF = [ RF ; A ];
        if feof(fid)
            fileEnd=true;
        end
        linenum = linenum + 1;
        line=fgets(fid);
        basisused_index=findstr(line,'$BASIS');
    end
    specs=RF(1:2:end) + 1i*RF(2:2:end);
    
    % GO 2022/01/24
    % LCModel uses negative BADELT values to encrypt basis sets
    % (LCModel.f source code lines 3666 and following)
    if dwelltime < 0
        dix = 1499;
        for rr = 1:length(specs)
            [randomresult, dix] = lcmodelrng(dix);
            specs(rr) = -specs(rr) .* exp(-20*randomresult + 10);
        end
    end
    
    specs=fftshift(conj(specs),1);
    vectorsize=length(specs);
    sz=[vectorsize 1];
    if mod(vectorsize,2)==0
        %disp('Length of vector is even.  Doing normal conversion');
        fids=fft(fftshift(specs,1),[],1);
    else
        %disp('Length of vector is odd.  Doing circshift by 1');
        fids=fft(circshift(fftshift(specs,1),1),[],1);
    end
    f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
    ppm=-f/(Bo*42.577);
    ppm=ppm+4.65;
    t=[dwelltime:dwelltime:vectorsize*dwelltime];
    txfrq=hzpppm*1e6;
    eval(['out.' metabName '.fids=fids;']);
    eval(['out.' metabName '.specs=specs;']);
    eval(['out.' metabName '.sz=[vectorsize 1 1 1];']);
    eval(['out.' metabName '.n=vectorsize;']);
    eval(['out.' metabName '.spectralwidth=spectralwidth;']);
    eval(['out.' metabName '.sz=sz;']);
    eval(['out.' metabName '.Bo=Bo;']);
    eval(['out.' metabName '.te=te;']);
    eval(['out.' metabName '.tr=[];']);
    eval(['out.' metabName '.dwelltime=1/spectralwidth;']);
    eval(['out.' metabName '.linewidth=linewidth;']);
    eval(['out.' metabName '.ppm=ppm;']);
    eval(['out.' metabName '.t=t;']);
    eval(['out.' metabName '.txfrq=txfrq;']);
    eval(['out.' metabName '.date=date;']);
    eval(['out.' metabName '.seq='''';']);
    eval(['out.' metabName '.sim='''';']);
    eval(['out.' metabName '.dims.t=1;']);
    eval(['out.' metabName '.dims.coils=0;']);
    eval(['out.' metabName '.dims.averages=0;']);
    eval(['out.' metabName '.dims.subSpecs=0;']);
    eval(['out.' metabName '.dims.extras=0;']);
    eval(['out.' metabName '.averages=1;']);
    eval(['out.' metabName '.flags.writtentostruct=1;']);
    eval(['out.' metabName '.flags.gotparams=1;']);
    eval(['out.' metabName '.flags.leftshifted=1;']);
    eval(['out.' metabName '.flags.filtered=0;']);
    eval(['out.' metabName '.flags.zeropadded=0;']);
    eval(['out.' metabName '.flags.freqcorrected=0;']);
    eval(['out.' metabName '.flags.phasecorrected=0;']);
    eval(['out.' metabName '.flags.averaged=1;']);
    eval(['out.' metabName '.flags.addedrcvrs=1;']);
    eval(['out.' metabName '.flags.subtracted=1;']);
    eval(['out.' metabName '.flags.writtentotext=1;']);
    eval(['out.' metabName '.flags.downsampled=0;']);
    eval(['out.' metabName '.flags.isISIS=0;']);
    
end

fclose(fid);
end
% RF=RF';
% rf(:,1)=RF(:,2)*180/pi;
% rf(:,2)=RF(:,1);
% rf(:,3)=ones(length(RF(:,1)),1);




% GO 2022/01/24
% LCModel uses negative BADELT values to encrypt basis sets
% (LCModel.f source code lines 3666 and following)
%++++++++++++++++ doubleprecision VERSION 2DP (MAR 1984) ++++++++++++++    4222
%  FUNCTION RANDOM.  PRODUCES A PSEUDORANDOM REAL ON THE OPEN INTERVAL      4223
%      (0.,1.).                                                             4224
%  DIX (IN doubleprecision) MUST BE INITIALIZED TO A WHOLE NUMBER          4225
%      BETWEEN 1.0D0 AND 2147483646.0D0 BEFORE THE FIRST CALL TO RANDOM       4226
%      AND NOT CHANGED BETWEEN SUCCESSIVE CALLS TO RANDOM.                  4227
%  BASED ON L. SCHRAGE, ACM TRANS. ON MATH. SOFTWARE 5, 132 (1979).         4228
%-----------------------------------------------------------------------    4229
function [randomresult,dix]=lcmodelrng(dix);

a   = 16807.0d0;
b15 = 32768.d0;
b16 = 65536.0d0;
p   = 2147483647.0d0;

 %                                                                           4231
 %  PORTABLE RANDOM NUMBER GENERATOR                                         4232
 %   USING THE RECURSION                                                     4233
 %    DIX = DIX*A MOD P                                                      4234
 %                                                                           4235
 %                                                                           4237
 %  7**5, 2**15, 2**16, 2**31-1                                              4238
                                                             
 %  GET 15 HI ORDER BITS OF DIX                                              4241
 xhi = dix./b16;
 xhi = xhi - rem(xhi,1.0d0);
 %  GET 16 LO BITS IF DIX AND FORM LO PRODUCT                                4244
 xalo =(dix-xhi.*b16).*a;
 %  GET 15 HI ORDER BITS OF LO PRODUCT                                       4246
 leftlo = xalo./b16;
 leftlo = leftlo - rem(leftlo,1.0d0);
 %  FORM THE 31 HIGHEST BITS OF FULL PRODUCT                                 4249
 fhi = xhi.*a + leftlo;
 %  GET OVERFLO PAST 31ST BIT OF FULL PRODUCT                                4251
 k = fix(fhi./b15);
 k = fix(k - rem(k,1.0d0));
 %  ASSEMBLE ALL THE PARTS AND PRESUBTRACT P                                 4254
 %   THE PARENTHESES ARE ESSENTIAL                                           4255
 dix =(((xalo-leftlo.*b16)-p)+(fhi-k.*b15).*b16) + k;
 %  ADD P BACK IN IF NECESSARY                                               4257
 if(dix < 0.0d0)
    dix = dix + p;
 end
 %  MULTIPLY BY 1/(2**31-1)                                                  4259
 randomresult = dix.*4.656612875d-10;
end %function random


