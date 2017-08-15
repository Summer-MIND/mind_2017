
clc
clear 
close all

%%

load('./files/MdPRLAnalysisExp3')
test                = str2func('signrank') ;

% for the subset of subjects with fewer number of estimations
probeTrialsAll      = [21 42 63 84 140 210 252 280, 301 322 343 364 420 490 532 560] ;
% initilaizing
rObj               = nan*ones(length(probeTrialsAll), length(subjects)) ;
rFt                = rObj ;

%%

for cnt_probe = 1:length(probeTrialsAll)
    pEstAll{cnt_probe}   = nan*ones(length(subjects),8) ;
end

clear sesdata
for cnt_sbj = 1:length(subjects)
    % load input and subjects' data
    inputname    = ['./inputs/input_',    subjects{cnt_sbj} , '.mat'] ; 
    resultsname  = ['./SubjectData/PRL_', subjects{cnt_sbj} , '.mat'] ;
    
    load(inputname)
    load(resultsname)
    
    probStim{1}     = expr.prob{1} ;
    probStim{2}     = expr.prob{2} ;
    cntprobe_half   = length(results.probEst)/2 ;
    cnt_probeMap    = [1:length(probeTrialsAll)] ;

    for cnt_probe = 1:length(results.probEst)
        
        % read the estimated Odds ratio and convert to probability
        probEst     = results.probEst{cnt_probe} ;
        probEst     = probEst./(1+probEst) ;
        cnt_probeAll= find(probeTrialsAll==expr.trialProbe(cnt_probe)) ;
        
        % read the actual probability 
        pObj         = probStim{2-mod(ceil(expr.trialProbe(cnt_probe)/expr.NtrialsShort),2)} ;
        pObj(pObj==0)= nan ;                             % remove the 9th object
        pObj         = round(pObj/0.05)*0.05 ;           % round up the numbers
        

        % for those who started with expr.flaginf==2
        if (expr.flaginf==1 && cnt_probe<=cntprobe_half) || (expr.flaginf==2 && cnt_probe>cntprobe_half)
        else
            probEst         = probEst' ;
            pObj            = pObj' ;
        end
        
        % for those that dont have all the estimation sessions
        if ~isempty(cnt_probeAll)
            pFt                                = myLL(ones(3,1)*nanmean(pObj,1), nanmean(pObj,2)*ones(1,3)) ;
            pEstAll{cnt_probeAll}(cnt_sbj,:)   = probEst(1:8) ;
            pObjAll{cnt_probeAll}(cnt_sbj,:)   = pObj(1:8) ;
            pFtAll{cnt_probeAll}(cnt_sbj,:)    = pFt(1:8) ;
            [rObj(cnt_probeAll, cnt_sbj), prObj(cnt_probeAll, cnt_sbj)]  = corr(probEst(1:8)' , pObj(1:8)') ;
            [rFt(cnt_probeAll, cnt_sbj),  prFt(cnt_probeAll, cnt_sbj)]   = corr(probEst(1:8)' , pFt(1:8)') ;
        end
    end
end

%%

% remove those that correlation is not significant
rObj(prObj>0.05 & prFt>0.05)    = nan ;
rFt(prObj>0.05 & prFt>0.05)     = nan ;
R                               = sum((rObj-rFt)'>0)./sum(~isnan(rObj')) ;
R                               = nanmean([R(1:8); R(9:end)])' ;

figure(1)
hold on
p = polyfit(probeTrialsAll(1:8)',R, 1) ;
plot(probeTrialsAll(1:8),R,'d', 'color', 'r','LineWidth',2)
set(gca,'FontName','Helvetica','FontSize',23,'FontWeight','normal','LineWidth',2,'yTick',0:0.2:1,'Xtick',[1 100:100:300], 'tickdir', 'out')
box off
axis([0 300 0.2 1])

% fit an exponential to data
x0                                     = [1*rand(1), 1*rand(1), 500*rand(1)];
lb                                     = [0, 0, 0]; 
ub                                     = [1, 1, 500];
[parE, resnorm]                        = lsqcurvefit(@ffitexp, x0,probeTrialsAll(1:8),R',lb,ub);
Xfit                                   = 1:280 ;
Yfit                                   = ffitexp(parE,Xfit) ;
plot(Xfit,Yfit,'--', 'color', 0.8*[1 0 0],'LineWidth',2) 
xlabel('trial (within a session)')
ylabel('fraction subjects')

cd ./figures
FigW = 6;
FigH = 5;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',[0,0,FigW,FigH],'units','centimeters');  
print('-dpdf','-r500','RindexExp3.pdf')
cd ../

%%

clear b
for cnt_probe = 1:length(probeTrialsAll)/2
    X            = [pEstAll{cnt_probe}; pEstAll{8+cnt_probe}] ;
    Y1           = [pObjAll{cnt_probe}; pObjAll{8+cnt_probe}] ;
    Y2           = [pFtAll{cnt_probe};  pFtAll{8+cnt_probe}] ;
    X            = X(:) ;
    Y            = [ones(size(Y1(:))) Y1(:) Y2(:)] ;

    b(:,cnt_probe)  = regress(X,Y) ;
    yCalc1          = Y*b(:,cnt_probe) ;
    y               = X ;
    Rsq(:,cnt_probe)= 1 - nansum((y - yCalc1).^2)/nansum((y - nanmean(y)).^2) ;
    res1(cnt_probe) = nansum(X-Y*b(:,cnt_probe)) ;
end

figure(2)
hold on
plot(probeTrialsAll(1:8), (b(2,:)./(sum(b(end-1:end,:))))','dr', 'LineWidth',2)
plot(probeTrialsAll(1:8), Rsq','dk', 'LineWidth',2)

% fit an exponential to data
x0                                     = [1*rand(1), 1*rand(1), 500*rand(1)];
lb                                     = [0, 0, 0]; 
ub                                     = [1, 1, 500];
[parE, resnorm]                        = lsqcurvefit(@ffitexp, x0,probeTrialsAll(1:8),Rsq,lb,ub);
Xfit                                   = 1:280 ;
Yfit                                   = ffitexp(parE,Xfit) ;
plot(Xfit,Yfit,'--', 'color', 0.5*[1 1 1],'LineWidth',2) 

x0                                     = [1*rand(1), 1*rand(1), 500*rand(1)];
lb                                     = [0, 0, 0]; 
ub                                     = [1, 1, 500];
[parE, resnorm]                        = lsqcurvefit(@ffitexp, x0,probeTrialsAll(1:8),(b(2,:)./(sum(b(end-1:end,:)))),lb,ub);
Xfit                                   = 1:280 ;
Yfit                                   = ffitexp(parE,Xfit) ;
plot(Xfit,Yfit,'--', 'color', 0.8*[1 0 0],'LineWidth',2) 

set(gca,'FontName','Helvetica','FontSize',23,'FontWeight','normal','LineWidth',2,'yTick',0:0.2:1,'Xtick',[1 100:100:300], 'tickdir', 'out')
box off
axis([0 300 0 1])
xlabel('trial (within a session)')
ylabel('relative weights, R^2')

cd ./figures
FigW = 6;
FigH = 5;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',[0,0,FigW,FigH],'units','centimeters');  
print('-dpdf','-r500','WindexExp3.pdf')
cd ../
