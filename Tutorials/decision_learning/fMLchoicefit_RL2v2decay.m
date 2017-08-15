function loglikehood = fMLchoicefit_RL2v2decay(xpar, sesdata)
% DESCRIPTION: fits data to object-based RL(2) model with decay using ML method
% INPUT: 
% sesdata structure which includes input, experiment and behavioral data
% OUTPUT:
% fitted parametres
% Version History:
% 0.1:  [2016-09-13]


% initializing parameters
loglikehood          = 0 ;
xpar([4:6])          = 1./(1+exp(-(xpar([4:6]))./sesdata.sig) ) ;      % make sure learning rates are between [0 1]
BiasL                = xpar(1) ;      % parameter to capture location bias 
magPatt              = xpar(2) ;      % subjective weight of pattern
magShape             = xpar(3) ;      % subjective weight of color
decay                = xpar(4) ;
alpha_rewPatt        = xpar(5) ;
alpha_rewShape       = xpar(5) ; 
alpha_unrPatt        = xpar(6) ;
alpha_unrShape       = xpar(6) ;

% reading session data
shapeMap             = sesdata.expr.targetShape ;
pattMap              = sesdata.expr.targetColor+3 ;
inputTarget          = sesdata.input.inputTarget ;
correcttrials        = sesdata.results.reward ;
choicetrials         = sesdata.results.choice ;
ntrials              = length(choicetrials) ;

% initilazing RL(2)
% v(1:3): value of the three shapes
% v(4:6): value of the three patterns
v                    = (0.5*ones(6,1)) ; 

% looping through trials
for cnt_trial=1:ntrials
    % after the break (mid experiment) informative and non-informative features change
    if cnt_trial==(1+(ntrials/2))
        magPatt     = xpar(3) ;
        magShape    = xpar(2) ;
    end
    
    correct         = correcttrials(cnt_trial) ;
    choice          = choicetrials(cnt_trial) ; 
    
    % identifying shape and pattern of the presented targets
    % 1: Left option, 2: Right option
    idx_shape(2)    = shapeMap(inputTarget(2, cnt_trial)) ;
    idx_patt(2)     =  pattMap(inputTarget(2, cnt_trial)) ;
    idx_shape(1)    = shapeMap(inputTarget(1, cnt_trial)) ;
    idx_patt(1)     =  pattMap(inputTarget(1, cnt_trial)) ;
    
    % estimating probability of choosing right option
    pChoiceR          = 1./(1+exp(-( magShape*(v(idx_shape(2))-v(idx_shape(1))) + ...
                                      magPatt*(v(idx_patt(2)) -v(idx_patt(1))) + BiasL ) )) ;
    pChoiceL          = 1-pChoiceR ;
    if choice == 2 
        loglikehood   = loglikehood - log(pChoiceR) ;
    else
        loglikehood   = loglikehood - log(pChoiceL) ; 
    end                      
    
    % updating estimates of value in RL(2)
    % updating the values of choosen target
    if correct
        idxC          = idx_patt(choice) ;
        v             = update(v, idxC, [], alpha_rewPatt) ;                   % potentiate the pattern selected 
        v             = decayV(v, 3+find([4:6]~=idx_patt(choice)), decay) ;    % decay the patterns not selected 

        idxC          = idx_shape(choice) ;
        v             = update(v, idxC, [], alpha_rewShape) ;                  % potentiate the shape selected 
        v             = decayV(v, find([1:3]~=idx_shape(choice)), decay) ;     % decay the shapes selected 
    else
        idxW          = idx_patt(choice) ;
        v             = update(v, [], idxW, alpha_unrPatt) ;                   % depress the pattern selected 
        v             = decayV(v, 3+find([4:6]~=idx_patt(choice)), decay) ;    % decay the patterns selected 

        idxW          = idx_shape(choice) ;
        v             = update(v, [], idxW, alpha_unrShape) ;                  % depress the shape selected 
        v             = decayV(v, find([1:3]~=idx_shape(choice)), decay) ;     % decay the shapes selected 
    end
    V(:,cnt_trial) = v ;
end % trial loop
end % funtion


function v = decayV(v, unCh, decay)
    v(unCh) = v(unCh) - (v(unCh)-0.5)*(decay) ;
end

function v = update(v, idxC, idxW, Q)
    if isempty(idxW) && ~isempty(idxC)
        v(idxC)   = v(idxC) + (1-v(idxC)).*Q ;
    elseif isempty(idxC) && ~isempty(idxW)
        v(idxW)   = v(idxW) - (v(idxW).*Q) ;
    elseif ~isempty(idxW) && ~isempty(idxC)
        v(idxC)   = v(idxC) + (1-v(idxC)).*Q ;
        v(idxW)   = v(idxW) - (v(idxW).*Q) ;
    elseif isempty(idxW) && isempty(idxC)
    end
end