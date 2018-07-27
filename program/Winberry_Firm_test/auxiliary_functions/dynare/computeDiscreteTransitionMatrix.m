function mTransition = computeDiscreteTransitionMatrix(vCapitalAdjust,vCapitalConstrained,vCutoff)
    
% Computes discrete transition matrix associated with decision rules following Young (2010)
% 
% Inputs
% (1) vCapitalAdjust: investment decision, conditional on adjusting (over fine grid)
%	(2) vCapitalConstrained: investment decision, conditional on not adjusting (over fine grid)
%	(3) vCutoff: fixed cost cutoff for adjusting capital (over fine grid)
%
% Outputs
% (1) mTransition: transition matrix
% 
% Thomas Winberry, Feburary 14th, 2018

% Declare global variables used in this function
global vCapitalGridFine nProdFine nCapitalFine nStateFine mProdTransition ppsiCapital
    
	
%----------------------------------------------------------------
% Compute weights for investment decision, conditional on adjusting
%----------------------------------------------------------------

% Compute nearest neighbor
vIndices 			= knnsearch(vCapitalGridFine,vCapitalAdjust);
vGridIndices 		= vCapitalGridFine(vIndices);

% Find indices and gridpoints above and below
vIndicesBelow 	= vIndices; vIndicesAbove = vIndices;

vIndicesBelow(vGridIndices > vCapitalAdjust) 	= vIndicesBelow(vGridIndices > vCapitalAdjust) - 1;
vIndicesBelow(vIndicesBelow < 1) 					= 1;
vGridBelow 													= vCapitalGridFine(vIndicesBelow);

vIndicesAbove(vGridIndices <= vCapitalAdjust) 	= vIndicesAbove(vGridIndices <= vCapitalAdjust) + 1;
vIndicesAbove(vIndicesAbove > nCapitalFine)		= nCapitalFine;
vGridAbove 													= vCapitalGridFine(vIndicesAbove);

% Compute weighting matrix
vAdjustWeightBelow 											= (vGridAbove - vCapitalAdjust) ./ (vGridAbove - vGridBelow);
vAdjustWeightBelow(vCapitalAdjust < vGridBelow) 	= 1;
vAdjustWeightBelow(vCapitalAdjust > vGridAbove) 	= 0;

vAdjustWeightAbove 											= (vCapitalAdjust - vGridBelow) ./ (vGridAbove - vGridBelow);
vAdjustWeightAbove(vCapitalAdjust < vGridBelow) 	= 0;
vAdjustWeightAbove(vCapitalAdjust > vGridAbove) 	= 1;

% Rename indices to denote conditional on adjusting
vIndicesAdjustBelow = vIndicesBelow; vIndicesAdjustAbove = vIndicesAbove;
clear vIndicesBelow vIndicesAbove


%----------------------------------------------------------------
% Compute weights for investment decision, conditional on not adjusting
%----------------------------------------------------------------

% Compute nearest neighbor
vIndices 		= knnsearch(vCapitalGridFine,vCapitalConstrained);
vGridIndices 	= vCapitalGridFine(vIndices);

% Find indices and gridpoints above and below
vIndicesBelow = vIndices; vIndicesAbove = vIndices;

vIndicesBelow(vGridIndices > vCapitalConstrained) 	= vIndicesBelow(vGridIndices > vCapitalConstrained) - 1;
vIndicesBelow(vIndicesBelow < 1) 						= 1;
vGridBelow 														= vCapitalGridFine(vIndicesBelow);

vIndicesAbove(vGridIndices <= vCapitalConstrained) 	= vIndicesAbove(vGridIndices <= vCapitalConstrained) + 1;
vIndicesAbove(vIndicesAbove > nCapitalFine) 			= nCapitalFine;
vGridAbove 															= vCapitalGridFine(vIndicesAbove);

% Compute weighting matrix
vConstrainedWeightBelow 													= (vGridAbove - vCapitalConstrained) ./ (vGridAbove - vGridBelow);
vConstrainedWeightBelow(vCapitalConstrained < vGridBelow) 	= 1;
vConstrainedWeightBelow(vCapitalConstrained > vGridAbove) 	= 0;

vConstrainedWeightAbove 													= (vCapitalConstrained - vGridBelow) ./ (vGridAbove - vGridBelow);
vConstrainedWeightAbove(vCapitalConstrained < vGridBelow) 	= 0;
vConstrainedWeightAbove(vCapitalConstrained > vGridAbove) 	= 1;

% Rename indices to denote conditional on not adjusting
vIndicesConstrainedBelow = vIndicesBelow; vIndicesConstrainedAbove = vIndicesAbove;
clear vIndicesBelow vIndicesAbove

%----------------------------------------------------------------
% Compute Transition Matrix
%----------------------------------------------------------------

% Compute matrices which record how many firms go from one capital
% state to another, splitting between gird points above and below
mAdjustTransitionBelow 			= zeros(nStateFine,nCapitalFine);
mAdjustTransitionAbove 			= zeros(nStateFine,nCapitalFine);
mConstrainedTransitionBelow 	= zeros(nStateFine,nCapitalFine);
mConstrainedTransitionAbove 	= zeros(nStateFine,nCapitalFine);
    
for k = 1:nCapitalFine

    mAdjustTransitionBelow(vIndicesAdjustBelow == k,k) =...
        vAdjustWeightBelow(vIndicesAdjustBelow == k);
    mAdjustTransitionAbove(vIndicesAdjustAbove == k,k) =...
        vAdjustWeightAbove(vIndicesAdjustAbove == k);
    mConstrainedTransitionBelow(vIndicesConstrainedBelow == k,k) =...
        vConstrainedWeightBelow(vIndicesConstrainedBelow == k);
    mConstrainedTransitionAbove(vIndicesConstrainedAbove == k,k) =...
        vConstrainedWeightAbove(vIndicesConstrainedAbove == k);
        
end

% Compute matrices which combine above and below, and extend along productivity dimension
mAdjustTransition 			= zeros(nStateFine,nStateFine);
mConstrainedTransition 	= zeros(nStateFine,nStateFine);

for iState = 1:nStateFine
    
   mAdjustTransition(iState,:) = reshape(ones(nProdFine,1) * (...
       mAdjustTransitionBelow(iState,:) + mAdjustTransitionAbove(iState,:)),...
       1,nStateFine);
   mConstrainedTransition(iState,:) = reshape(ones(nProdFine,1) * (...
        mConstrainedTransitionBelow(iState,:) + mConstrainedTransitionAbove(iState,:))...
        ,1,nStateFine);
    
end

% Combine the adjust and constrained decision, and weight by productivity draws
mProdTransitionExpanded 	= repmat(mProdTransition,nCapitalFine);
mCutoff 							= vCutoff * ones(1,nStateFine);
mTransition 						= ((mCutoff ./ ppsiCapital) .* mAdjustTransition +...
											(1 - (mCutoff ./ ppsiCapital)) .* mConstrainedTransition) .* ...
											mProdTransitionExpanded;
