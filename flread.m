clear all; close all;
[Y1,FS1,NBITS1,OPTS1]=wavread('sgbb1.wav');
x1=Y1;
% [P,Q] = rat(44100/FS1);
% x1 = resample(Y1,P,Q);
% left = Y1(:,1);
% right = Y1(:,2);
% x1 = 0.5*(left+right);
% FS1=44100;

% player = audioplayer(x1,FS1);
% play(player);
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

%%%fid = fopen('sgbb1_us.txt', 'a+');
th=0;
for i=1:numOfFrames-1
    xwn = (x1(curPos:curPos+floor(wintime*FS1)-1));
    [px,f] = pwelch(xwn,FS1);
    %maxfreq=freq(find(pxx>threshold))
    th=ceil(max(f-px))/10;
end
th
for i=1:numOfFrames-1
%     disp(i);
    cnt=0;
    cntneg = 0;
    cntpos = 0;
    xwin = (x1(curPos:curPos+floor(wintime*FS1)-1));
%     c=abs(rhythm(xwin));
%     figure; plot(c)
    startTime=curPos/FS1;
%     length(smooth(c(:,4)))
    endTime=(curPos+floor(wintime*FS1)-1)/FS1;
    
    
% % %     pxx = periodogram(xwin);
    [pxx,f] = pwelch(xwin,FS1);
    %maxfreq=freq(find(pxx>threshold))
    %th=ceil(max(f-pxx))/10
    
    

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
%           slope(:,i) = (x1(j) -x1(j-1))
            if ((c1(j)-c1(j-1)) == 0)
                cnt=cnt+1;
            end
            if ((c1(j)-c1(j-1)) < 0)
                cntneg=cntneg+1;
            
            end
            if ((c1(j)-c1(j-1)) > 0)
                cntpos=cntpos+1;
                muposfrms=[muposfrms mean(xwin)];
%               dat=[dat; xwin];
            end
        end
        if (cnt >0)
            frms=[frms cnt];
%         cnt
    
        end   
        mu=mean(xwin);
        sigma=std(xwin);
        
        if (cntneg >0)
            frmsneg=[frmsneg cntneg];
            munegfrms=[munegfrms mu];
            signegfrms=[signegfrms sigma];
%         cntneg
%         figure; plot(smooth(c(:,4)))
        end
    
        if (cntpos >0)
            frmspos=[frmspos cntpos];
%         cntpos
%         figure; plot(smooth(c(:,4)))
        end
        timeInterval1(NZcntr) = startTime;
        timeInterval2(NZcntr) = endTime;
%         fprintf(fid, '%d\t%.3f\t%.3f\n', i, startTime, endTime);        
    end
%     [ymax,imax,ymin,imin] = extrema2(smooth(c(:,4)));
%     hold on;
%     plot(smooth(c(imax)),ymax,'r*',smooth(c(imin)),ymin,'g*')
%     [pks,locs] = findpeaks(smooth(c(:,4)))
%     upperlimit = ones(1, 1)*1;
%     plot(upperlimit);
%     hold on;
%     plot(upperlimit*0);
%     for j=1:6
%         c1 =abs(c(:,3));
%         c1 = c1/max(c1);
        
%     end
    
    
%     for j =1:6
% %         c = c/max(c);
%         figure; plot(c)
%     end
%     pause(0.4);
    curPos = curPos + floor(steptime*FS1);
end
% fclose(fid);
[s,indx] = sort(frms);
[s1,indx1] = sort(frmsneg);
[s2,indx2] = sort(frmspos);

[smu1,idx1] = sort(signegfrms);
munegfrms
signegfrms

if (length(signegfrms)>0)

    [h,p,ci,stats] = ttest(signegfrms);

% n=kmeans(s1',4)
% 
% gmfit = gmdistribution.fit(s1',4);
% 
% [idx,nlogl] = cluster(gmfit,s1')

% d = pdist(smu1');
% d = pdist(smu1');
% Z = linkage(d,'single')
% cl = cluster(Z,'maxclust',2:8)
% 
% dendrogram(Z)
% 
% C=nancov(s1');
% D = pdist(s1','mahalanobis',C)
% 
% [kclidx2, kclval2]= kmeans(s1', 2);
% [kclidx3, kclval3]= kmeans(s1', 3);
% [kclidx4, kclval4]= kmeans(s1', 4);
% [kclidx5, kclval5]= kmeans(s1', 5);
% [kclidx6, kclval6]= kmeans(s1', 6);
% [kclidx7, kclval7]= kmeans(s1', 7);
% 
% C=nancov(smu1');
% D = pdist(smu1','mahalanobis',C);
% 
% [kclidx2, kclval2]= kmeans(smu1', 2);
% [kclidx3, kclval3]= kmeans(smu1', 3);
% [kclidx4, kclval4]= kmeans(smu1', 4);
% [kclidx5, kclval5]= kmeans(smu1', 5);
% [kclidx6, kclval6]= kmeans(smu1', 6);
% [kclidx7, kclval7]= kmeans(smu1', 7);
% 
% kcl(1,:)=kclidx2;
% kcl(2,:)=kclidx3;
% kcl(3,:)=kclidx4;
% kcl(4,:)=kclidx5;
% kcl(5,:)=kclidx6;
% kcl(6,:)=kclidx7;

% kcl'

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

% for i=1:size(cl,1)
%     mode(cl(i,:))
% end

% w=[1 1 1];
% votes=cl;
% options=['1','2','3','4']';
% OPTIONS=repmat(options,[1,size(w,2),size(votes,1)]);
% B=bsxfun(@eq, permute(votes,[3 2 1]),OPTIONS);
% W=squeeze(sum(bsxfun(@times,repmat(w,size(options,1),1),B),2))';
% [xx,i]=max(W,[],2);
% options(i)

% [d,n] = size(cl(1,:)');
% w=ones(1,n);

% cl(1,:)'

% E1= floor(sum(cl(1,:)')/length(cl(1,:)));
% E2= floor(sum(cl(2,:)')/length(cl(2,:)));
% E3= floor(sum(cl(3,:)')/length(cl(3,:)));
% E4= floor(sum(cl(4,:)')/length(cl(4,:)));
% E5= floor(sum(cl(5,:)')/length(cl(5,:)));
% E6= floor(sum(cl(6,:)')/length(cl(6,:)));
% E7= floor(sum(cl(7,:)')/length(cl(7,:)));
% E1
% E2
% E3
% E4
% E5
% E6
% E7

% % function [value,vote] = majority(t1,w);
% % % MAJORITY : returns (weighted) majority vote
% % % [value,vote] = majority(t1[,w])
% % %	t1    - the data
% % %                if vector, returns scalar
% % %                if matrix, returns row vector (i.e. majority over columns)
% % %	w     - the weight vector, if omitted equal weights
% % %	value - the result (e.g. majority([1 1 2 3]) = 1)
% % %	vote  - the vote supporting the result (for above example, vote = 2)
% % 
% % % Copyright (c) 1995 Frank Dellaert
% % % All rights Reserved
% % 
% % [d,n] = size(t1);
% % if (nargin==1), w=ones(1,n); end
% % 
% c2=cl(7,:);
% if(d==1)
%   [value,vote] = majority1(cl',w);
% else
%   w=ones(1,d); 
%   value = zeros(1,n);
%   vote = zeros(1,n);
%   for i=1:n
%     [value(i),vote(i)] = majority1(cl(:,i),w)
%   end
% end
% vote
% value

% save('data','x1','-ascii')
% load('data')
% label = emgm(data',4); 
% spread(data',label);

xwin = x1(curPos:L-1);
% ES(:,numOfFrames) = rhythm(xwin);
% rhythm(xwin)
