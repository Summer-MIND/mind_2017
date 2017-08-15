
clc
clear 
close all

%%

load('./files/MdPRLAnalysisExp3') ;

RL2_decay           = cat(1,mlparRL2_decay{:}) ;
RL2conj_decay       = cat(1,mlparRL2conj_decay{:}) ; 
test                = str2func('signrank') ;

%%

clear sesdata
for cnt_sbj = 1:length(subjects)
    % load input and subjects' data
    inputname    = ['./inputs/input_',    subjects{cnt_sbj} , '.mat'] ; 
    resultsname  = ['./SubjectData/PRL_', subjects{cnt_sbj} , '.mat'] ;
    
    load(inputname)
    load(resultsname)

    % generate fitting session data
    sesdata.sig                   = 0.2 ;
    sesdata.input                 = input ;
    sesdata.expr                  = expr ;
    sesdata.results               = results ;

    % fitting model: RL2 decay
    xpar                             = RL2_decay(cnt_sbj, 1:6) ;
    LLAll_RL2(cnt_sbj,:)             = fMLchoiceLL_RL2v2decay(xpar, sesdata) ;
    
    % RL2 conjunction decay
    xpar                             = RL2conj_decay(cnt_sbj, 1:5) ;
    LLAll_RL2conj(cnt_sbj,:)         = fMLchoiceLL_RL2conjdecay(xpar, sesdata) ;
end

%%

wSize               = 20 ;
% average each subject's first half and second half
LL1                 = cat(2, LLAll_RL2') ;
LL1                 = mean(reshape(LL1, [], 2, length(subjects)),2) ;
LL1                 = reshape(LL1,[], length(subjects)) ;

LL2                 = cat(2, LLAll_RL2conj') ;
LL2                 = mean(reshape(LL2, [], 2, length(subjects)),2) ;
LL2                 = reshape(LL2,[], length(subjects)) ;

% filter data using a moving avreage box
uX                  = wSize:280 ;
X1                  = filter(ones(1,wSize)/wSize,1,LL1) ;
X1                  = X1(wSize:end,:) ;
mu1                 = mean(X1,2)' ;
sd1                 = std(X1')./sqrt(length(subjects)) ;

X2                  = filter(ones(1,wSize)/wSize,1,LL2) ;
X2                  = X2(wSize:end,:) ;
mu2                 = mean(X2,2)' ;
sd2                 = std(X2')./sqrt(length(subjects)) ;

% prepare for patch
x1                  = [uX fliplr(uX)];
y1                  = [mu1+sd1 fliplr(mu1-sd1)];

x2                  = [uX fliplr(uX)];
y2                  = [mu2+sd2 fliplr(mu2-sd2)];

figure(1)
hold on
hpatch    = patch(x1,y1,'b'); 
set(hpatch,'EdgeColor','none'); 
set(hpatch,'FaceColor',[0 .0 1]); 
hline     = plot(uX,mu1,'b-'); 
set(hline,'LineWidth',2); 
set(hline,'Color','b'); 
box off
alpha(hpatch,0.3);

hpatch    = patch(x2,y2,'r'); 
set(hpatch,'EdgeColor','none');
set(hpatch,'FaceColor',[1.0 .0 0]); 
hline     = plot(uX,mu2,'r-'); 
set(hline,'LineWidth',2); 
set(hline,'Color','r'); 
box off
alpha(hpatch,0.3);

plot(1:280, -log(0.5)*ones(1,280), '--k')
set(gca,'FontName','Helvetica','FontSize',25,'FontWeight','normal','LineWidth',2,'XTick',[1 100:100:350],...
        'ytick',0.0:0.1:0.7)
set(gca,'TickDir','out')
ylabel('goodness-of-fit')
xlabel('trial (within a session)')
axis([1 300 0.39 .7])

cd ./figures
FigW = 6;
FigH = 5 ;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',[0,0,FigW,FigH],'units','centimeters');  
print('-dpdf','-r500','LLtimeExp3.pdf')
cd ../

