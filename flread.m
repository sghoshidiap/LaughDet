%@COPYRIGHT Sucheta Ghosh sucheta.ghosh@idiap.ch
clear all; close all;
[Y1,FS1,NBITS1,OPTS1]=wavread('sgbb1.wav');
x1=Y1;

wintime = 2.5;
steptime = 1.25;
L = length(x1)
totdur = L/FS1;
numOfFrames = floor(totdur/steptime);
curPos = 1;
winpt = floor(wintime*FS1)
ES = zeros(floor(wintime*FS1),numOfFrames);
N1 = 0;
slope=[];

frms=[];
frmsneg=[];
frmspos=[];

mufrms=[];
munegfrms=[];
muposfrms=[];
dat=[];

signegfrms=[];
sigposfrms=[];

timeInterval1 = [];
timeInterval2 = [];
NZcntr=0;

th=0;
for i=1:numOfFrames-1
    xwn = (x1(curPos:curPos+floor(wintime*FS1)-1));
    [px,f] = pwelch(xwn,FS1);
    %maxfreq=freq(find(pxx>threshold))
    th=ceil(max(f-px))/10;
end
th
for i=1:numOfFrames-1
    cnt=0;
    cntneg = 0;
    cntpos = 0;
    xwin = (x1(curPos:curPos+floor(wintime*FS1)-1));
    startTime=curPos/FS1;
    endTime=(curPos+floor(wintime*FS1)-1)/FS1;
    
    [pxx,f] = pwelch(xwin,FS1);
    
    ytmp = pxx;
    ytmp(find(pxx>=th))= i;
    ytmp(find(pxx < th)) = 0;
    NZ = any(ytmp);
    mu=0; sigma=0;
    
    if (NZ)
        NZcntr=NZcntr+1;
        c=abs(rhythm(xwin));
        figure; plot(c(:,3))
        sprintf('%d\t%.3f\t%.3f\n', i, startTime, endTime);
        c1=c(:,3);
        for j = 2:length(c1)
            if ((c1(j)-c1(j-1)) == 0)
                cnt=cnt+1;
            end
            if ((c1(j)-c1(j-1)) < 0)
                cntneg=cntneg+1;
            
            end
            if ((c1(j)-c1(j-1)) > 0)
                cntpos=cntpos+1;
                muposfrms=[muposfrms mean(xwin)];
            end
        end
        if (cnt >0)
            frms=[frms cnt];
    
        end   
        mu=mean(xwin);
        sigma=std(xwin);
        
        if (cntneg >0)
            frmsneg=[frmsneg cntneg];
            munegfrms=[munegfrms mu];
            signegfrms=[signegfrms sigma];
        end
    
        if (cntpos >0)
            frmspos=[frmspos cntpos];
        end
        timeInterval1(NZcntr) = startTime;
        timeInterval2(NZcntr) = endTime;
    end

    curPos = curPos + floor(steptime*FS1);
end
[s,indx] = sort(frms);
[s1,indx1] = sort(frmsneg);
[s2,indx2] = sort(frmspos);

[smu1,idx1] = sort(signegfrms);
munegfrms
signegfrms

if (length(signegfrms)>0)

    [h,p,ci,stats] = ttest(signegfrms);


    threshold=ci(2)-stats.sd;
%     threshold=ci(2)

    for i=1:size(signegfrms,2)
        if (signegfrms(i)>threshold)
            i
            timeInterval1(i)
            timeInterval2(i)
            fprintf(fid, '%d\t%.3f\t%.3f\n', i, timeInterval1(i), timeInterval2(i));
        end
    end

end

fclose(fid);

xwin = x1(curPos:L-1);
