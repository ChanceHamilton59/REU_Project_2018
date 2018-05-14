function [thisStimImage, flipStateA, flipStateB, rSeedA, rSeedB] = stimQuilter(stim, orintation)
%STIMQUILTER Summary of this function goes here
%   Detailed explanation goes here

global shareMask
global carrAPattsGlobal 
global carrBPattsGlobal

env.type = 'halfDisc';
env.taper = 0.1;
env.lambda = [];
env.phase = 0;
env.shift = 0;
if orintation == 1
   env.ori = -45; 
else
   env.ori = 45; 
end


scr.gray = 0.5;
scr.inc  = 0.5;
scr.white = 1;
scr.black = 0;

carrPattA   = stim.carrPattA - mean(mean(stim.carrPattA));
carrPattB   = stim.carrPattB - mean(mean(stim.carrPattB)); 


imgSize     = size(carrPattA,1); 
envPatt     = stimMakeEnvPattPsy(env,imgSize); % range -1 to +1

winMask     = stimMakeCosTaper(imgSize,0.2);

flipStateA  = [];  
flipStateB  = [];   % (something to return, if these are not used)
rSeedA      = 0;       
rSeedB      = 0;





if strcmp(stim.modType,'quilt') ...
       || strcmp(stim.modType,'quilt_con') % quilts a la Landy-Oruc
    
    envPattA = sqrt((1 + (stim.modDepth/100)*envPatt)/2.0).*winMask;  % range 0 to +1
    envPattB = sqrt((1 - (stim.modDepth/100)*envPatt)/2.0).*winMask;    
    
    carrPattA = carrPattA./max(abs(carrPattA(:))); % range -1 to +1, with mean = 0.0 
    %carrPattA = (carrPattA + 1.0)/2;               % range 0 to +1, with mean = 0.5
    
    
    halfPattA = carrPattA.*envPattA;  % re-envelope, now with zero-mean -> range -.5 to +.5
    halfPattA = 2*halfPattA;
    
    
    carrPattB = carrPattB./max(abs(carrPattB(:))); % range -1 to +1, with mean = 0.0
    %carrPattB = (carrPattB + 1.0)/2.;               % range 0 to +1, with mean = 0.5
    

    
    % finally, join A to B:
    halfPattB = carrPattB.*envPattB;  % re-envelope, now with zero-mean -> range -.5 to +.5
    halfPattB = 2*halfPattB;
    
   
    
    thisStimImage =  halfPattA + halfPattB;  % quilt, range -1 to +1
       


    
    % thisStimImage = scr.gray + scr.inc*thisStimImage;
    
    %figure(99); imagesc(thisStimImage); axis square off; colormap('gray');
    
    
    




end

