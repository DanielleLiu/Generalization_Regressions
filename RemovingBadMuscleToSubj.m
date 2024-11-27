function normalizedTMFullAbrupt= RemovingBadMuscleToSubj(normalizedTMFullAbrupt, subjToRemove, muscleToRemove)
%Code to make sure we remove the same muscle for all the analysis


% Remmoved muscle removed the entire muscle for the data. This muscle are
% being removed due to bad normalization. 
%You can see the subjects and the correspondant muscle that is being remove

 % Treadmill group (Experiment 1)
[normalizedTMFullAbrupt]=RemoveBadMuscles(normalizedTMFullAbrupt,{'BATR13','BATR03'},{{'fVLs','sVLs','fVMs','sVMs'},{'fVLs','sVLs','fVMs','sVMs'}});

% Generalization group (Experiment 2)
 [normalizedTMFullAbrupt]=RemoveBadMuscles(normalizedTMFullAbrupt,{'BATS02','BATS04','BATS06','BATS09','BATS12'},...
     {{'fSOLs','sSOLs','fVMs','sVMs','fVLs','sVLs','sRFs','fRFs'},{'fBFs','sBFs'},{'fRFs','sRFs'},{'fRFs','sRFs'},{'sRFs','fRFs','fVLs','sVLs'}});
 
 if nargin > 2
     %If a list of muscle is provided.
     [normalizedTMFullAbrupt]=RemoveBadMuscles(normalizedTMFullAbrupt,subjToRemove,muscleToRemove);
 end
end
 

