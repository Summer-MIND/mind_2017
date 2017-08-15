
clc
clear 
close all

%%
% initializing the values for fit
nrep              = 5 ;                         % how many times to repeat each fit
op                = optimset;                   % generate optimization options structure
warning('off')
subjects          = {'ad01', 'ag01', 'ah01', 'al01', 'av01', 'aw01', 'bc01', 'bd01', ...
                     'bf01', 'bh01', 'bj01', 'bl01', 'bn01'} ;

%% loop over all the subjects
for cnt_sbj = 1:length(subjects)
    % load input and subjects' data
    inputname    = ['./inputs/input_',    subjects{cnt_sbj} , '.mat'] ; 
    resultsname  = ['./SubjectData/PRL_', subjects{cnt_sbj} , '.mat'] ;
    
    load(inputname)
    load(resultsname)
    
    flaginf(cnt_sbj,1)            = expr.flaginf ;                 % which dimension was more informative at first
    rew{cnt_sbj,1}                = results.reward ;

    % generate fitting session data
    sesdata.sig                   = 0.2 ;
    sesdata.input                 = input ;
    sesdata.expr                  = expr ;
    sesdata.results               = results ;
    
    % initialize fit performance
    fvalminRL2_decay              = length(sesdata.results.reward) ;
    fvalminRL2conj_decay          = length(sesdata.results.reward) ;

    % reapeat each fit to avoid finding local minima
    for cnt_rep = 1:nrep
        disp(['-----------------------------------------------------------'])
        disp(['Subject: ', num2str(cnt_sbj),', Repeat: ', num2str(cnt_rep)])

        %% fitting model: RL2 decay
        Nparam                           = 6 ;
        ipar                             = rand(1,Nparam)  ;
        [xpar fval exitflag output]      = fminsearch(@fMLchoicefit_RL2v2decay, ipar, op, sesdata) ;
        if fval <= fvalminRL2_decay
            xpar([4:6])                  = 1./(1+exp(-(xpar([4:6]))./sesdata.sig) ) ;
            fvalminRL2_decay             = fval ;
            mlparRL2_decay{cnt_sbj}(1:6) = xpar(1:6) ;
            mlparRL2_decay{cnt_sbj}(100) = fval ;
            mlparRL2_decay{cnt_sbj}(101) = fval./length(sesdata.results.reward) ;
            mlparRL2_decay{cnt_sbj}(102) = output.iterations;
            mlparRL2_decay{cnt_sbj}(103) = exitflag ;
        end
        
        %% fitting model: RL2 conjunction decay
        NparamBasic                      = 5 ;
        ipar                             = rand(1,Nparam)  ;
        [xpar fval exitflag output]      = fminsearch(@fMLchoicefit_RL2conjdecay, ipar, op, sesdata) ;
        if fval <= fvalminRL2conj_decay
            xpar([3:5])                      = 1./(1+exp(-(xpar([3:5]))./sesdata.sig) ) ;
            fvalminRL2conj_decay             = fval ;
            mlparRL2conj_decay{cnt_sbj}(1:5) = (xpar(1:5)) ;
            mlparRL2conj_decay{cnt_sbj}(100) = fval ;
            mlparRL2conj_decay{cnt_sbj}(101) = fval./length(sesdata.results.reward) ;
            mlparRL2conj_decay{cnt_sbj}(102) = output.iterations;
            mlparRL2conj_decay{cnt_sbj}(103) = exitflag ;
        end
    end
end
save ./files/MdPRLAnalysisExp3
