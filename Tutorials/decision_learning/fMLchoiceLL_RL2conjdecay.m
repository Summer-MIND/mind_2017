function loglikehood = fMLchoiceLL_RL2conjdecay(xpar, sesdata)
% DESCRIPTION: fits data to object-based RL(2) model with decay using ML method
% INPUT: 
% sesdata structure which includes input, experiment and behavioral data
% OUTPUT:
% fitted parametres
% Version History:
% 0.1:  [2016-09-13]


% initializing parameters
BiasL                = xpar(1) ;      % parameter to capture location bias 
mag                  = xpar(2) ;      % subjective weight
decay                = xpar(3) ;
alpha_rew            = xpar(4) ;
alpha_unr            = xpar(5) ; 

inputTarget          = sesdata.input.inputTarget ;
correcttrials        = sesdata.results.reward ;
choicetrials         = sesdata.results.choice ;
ntrials              = length(choicetrials) ;

% initilazing RL(2)
v                    = (0.5*ones(9,1)) ; 

for cnt_trial=1:ntrials
    
    correct           = correcttrials(cnt_trial) ;
    choice            = choicetrials(cnt_trial) ; 
    
    % estimating probability of choosing right option
    % 1: Left option, 2: Right option
    pChoiceR          = 1./(1+exp(-( mag*(v(inputTarget(2, cnt_trial))-v(inputTarget(1, cnt_trial))) + BiasL ) )) ;
    pChoiceL          = 1-pChoiceR ;
    if choice == 2 
        loglikehood(cnt_trial)   = - log(pChoiceR) ;
    else
        loglikehood(cnt_trial)   = - log(pChoiceL) ; 
    end                      
    
    % updating estimates of value in RL(2)
    if correct
        idxC          = inputTarget(choice, cnt_trial) ;
        v             = update(v, idxC, [], alpha_rew) ;                                   % potentiate the object selected 
        v             = decayV(v, find([1:9]~=inputTarget(choice, cnt_trial)), decay) ;    % decay the objects not selected 
    else
        idxW          = inputTarget(choice, cnt_trial) ;
        v             = update(v, [], idxW, alpha_unr) ;                                   % potentiate the object selected 
        v             = decayV(v, find([1:9]~=inputTarget(choice, cnt_trial)), decay) ;    % decay the objects not selected 
    end
    
    V(:,cnt_trial) = v ;
    
end
end

function v = decayV(v, unCh, decay)
    v(unCh) = v(unCh) - (v(unCh)-0.5)*(decay) ;
end

function v = update(v, idxC, idxW, Q)
    if isempty(idxW)
        v(idxC) = v(idxC) + (1-v(idxC)).*Q ;
    elseif isempty(idxC)
        v(idxW) = v(idxW) - (v(idxW).*Q) ;
    elseif ~isempty(idxW) && ~isempty(idxC)
        v(idxC) = v(idxC) + (1-v(idxC)).*Q ;
        v(idxW) = v(idxW) - (v(idxW).*Q) ;
    end
end
