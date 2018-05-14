function [thisStimImage flipStateA flipStateB rSeedA rSeedB bubbleMask] = stimMakeTexSegCD(stim,scr,winMask,diagnostics,preComp,pattInd)
%
% stimMakeTexSeg - create a stimulus image, for texture segregation experiment
%
% usage:  [thisStimImage flipStateA flipStateB rSeedA rSeedB] = stimMakeTexSeg(stim,scr,winMask,task.diag)


% TO DO:
%   diagnostic images, printouts:  make contingent on  task.diag.images == 1
%   incorporate option for stim.carrA/B.type = string (e.g. "pos_gaussian", "gabor") rather than integer
%       removal of carrier means contingent of "pos_gaussian" or "neg_gaussian" rather than 10 or 11
%   OR, incorporate another element in stim.carrA/B, for whether zero-balanced or not
%  =>in either case, put the mean-removal within stimMakeCarrPattPsy, contingent on 'synthTex' 



% questions:
%   image-quilts:  is clipping problem greater at smaller modDepths ? 
%   conMod or quilt:  still works ok if .type is NOT 'image' ?


% needs to run:
%   stimMakeEnvPattPsy.m - create envelope
%   stimMakeCarrPattPsy.m - create carrier (texture)
%       stimMake2dFractNoise.m
%           stimMakeFractSclMatrix.m
%       drawSynthTexReg.m
%       drawSynthTexNat.m
%       psyTestQuadStats.m
%       phaseScramble.m
%       imwhiten.m
%       imgZclip.m - clip extreme values
%   imgPrintStats.m


% related files:
%    stimMakeTexSeg_OLD.m - using "old" method of handling mean/contrast


% recent chagnes:
% 16-Nov-07, LA: 
% 	- added modulation type lumcon
% 	- debugged square frame in contrast modulation
% 	- added luminance profile to output figure if stim.modType ='lumcon'. 
% 17-21 Nov 07, CB, modType='con':
% 	- elimination of lum-artifact for modType='con
%   - incorporated added luminance-mod into modType='con'
%   - approx calibrated rms-contrast for modType='con'
% Jan 08 ?? - AY - new, improved way to control contrast/luminance
% 23-25 Apr 08, CB
%   - restored lumConAdd
%   - return flipState's, rSeed's provided by stimMakeCarrPattPsy
% (2008-2013, LA - modifications for thesis ...)
% 26 Sep 13, CB: ## prevent removal of carrier means for mpType==10 or 11
%       (not working yet, maybe incorrect !!)

global shareMask
global carrAPattsGlobal 
global carrBPattsGlobal

shareMask = winMask;

imgSize = stim.size;

bubbleMask = zeros(imgSize);

envPatt = stimMakeEnvPattPsy(stim.env,imgSize); % range -1 to +1

flipStateA=[];  flipStateB=[];   % (something to return, if these are not used)
rSeedA=0;       rSeedB=0;
% 
% bits.black = 0;  
% bits.white = 2^14; 
% bits.gray = (bits.black + bits.white) / 2; 
% bits.inc = abs(bits.white - bits.gray);

% if strcmp(stim.modType,'lum')   % luminance-modulation
%     envPatt = winMask.*envPatt;     % apply window mask     
%     thisStimImage = scr.gray + (stim.modDepth/100)*scr.inc *envPatt; % range 0 to 255          
% 
%     
% elseif strcmp(stim.modType,'con')   % contrast-modulation of carrier-A
%     if isfield(stim.env,'bubble');
%         envPatt = (envPatt + 1)./2; %range 0 to 1
%         envMask = ones(imgSize);
%         
%         bs = stim.env.bubble.s;
%         
%         % make bubble with gauss to size (env.bubble.s = diameter)
%         bubble = 1 - stimMakeGausWindow(bs,bs/5);
%     
%         % scatter env.bubble.n bubbles over the envelope
%         for i = 1:stim.env.bubble.n
%             posx = round(1 + (imgSize-bs)*rand(1));%randomize position
%             posy = round(1 + (imgSize-bs)*rand(1));
%               
%         envMask(posx:posx+bs-1, posy:posy+bs-1) = ...  %place micropattern in position   
%             envMask(posx:posx+bs-1, posy:posy+bs-1).*bubble;
%         end
%         
%         envMask = abs(envMask-1);
%         %envPatt =  envPatt .* envMask;
%         %envPatt = abs(envPatt-1);
%         
%         bubbleMask = envMask;
%         %bubbleMask = envPatt .* winMask; 
%     end
% 
%     envPatt = ((stim.modDepth/100)*envPatt + 1.0)/2.0;  % range 0 to 1
%     envPatt = envPatt.*winMask;     % apply mask now, since we care only about what will be on the screen
%        
%     if(preComp)
%         
%         if(isfield(stim,'AorB' ))
%         
%             if(stim.AorB)
%                 carrPattA  = double(carrAPattsGlobal{pattInd});
%                 carrPattB  = carrPattA;
%             else
%                 carrPattA  = double(carrBPattsGlobal{pattInd});
%                 carrPattB  = carrPattA;
%             end
%         
%         else
%             carrPattA  = double(carrAPattsGlobal{pattInd});
%             carrPattB  = carrPattA; 
%         end
%         
%         %    carrPattB  = double(carrBPattsGlobal{pattInd});
%         flipStateA = 'nn';
%     %    flipStateB = 'nn';
%         rSeedA     = ClockRandSeed;
%     %    rSeedB     = ClockRandSeed;        
%     else
%     
%     
%         [carrPattA flipStateA rSeedA] = stimMakeCarrPattPsy(stim.carrA,scr.gray,stim.size); 
%     % sdClipped, range = 0 to zMax;  for synthTex, range = 0-255
%         
%     end 
%         
% %     figure(2);  
% %     subplot (2,3,1); imagesc(carrPattA,[0 256]);colormap(gray);title('1st raw carrPattA'); colorbar;%$$-diag
% %     figure(3);
% %     subplot(2,3,1); plot(1:480,carrPattA(240,:), 1:480,127.5*ones(1,480));   % $$ diag  
%     
%     if ~(stim.carrA.mpType==10 || stim.carrA.mpType==11)   % ##-CB  i.e., if NOT a gaussian
%         carrPattA = carrPattA - mean(carrPattA(:));     % zero-mean
%     elseif (stim.carrA.mpType==10 || stim.carrA.mpType==11)  % if pos_gaussian or neg_gaussian
%         carrPattA = carrPattA - 127.5;  % $$ correct ?   if good, should use scr.gray ?
%     end
%     
%     % range now +/-127.5
% %     figure(2);
% %     subplot(2,3,2);imagesc(carrPattA,[-128 128]);colormap(gray);title('later raw carrPattA');colorbar;%$$ diag
% %     figure(3);
% %     subplot(2,3,2); plot(1:480,carrPattA(240,:), 1:480,0*ones(1,480));  % $$ diag
% 
%     if isfield(stim.env,'bubble');
%         carrPattA = carrPattA .* bubbleMask;
%         bubbleMask = bubbleMask .* winMask;
%     end
%     
%     carrPattA = carrPattA./max(abs(carrPattA(:))); % -> range -1 to +1, with mean = 0.0
%     carrPattA = (carrPattA + 1.0)./2;               % range 0 to +1, with mean = 0.5
%     
% %     figure(2);
% %     subplot(2,3,3);  imagesc(carrPattA,[0 1]); colormap(gray); title('adjusted carrPattA'); colorbar; %$$ diag
% %     figure(3);
% %     subplot(2,3,3); plot(1:480,carrPattA(240,:), 1:480,0.5*ones(1,480));    % $$ diag
% 
%     if ~strcmp(stim.carrA.type, 'synthNoOverlap')
%         cv = std2(carrPattA)/mean2(carrPattA);  % i.e., coefficient of variation, typically 0-1 
%                                         % (note denominator is ALWAYS = 0.5, due to above code)
%         conScl = stim.carrA.con / (100*cv);    % scale ot specified carrier contrast
%     
%         carrPattA = 0.5 + conScl*(carrPattA - 0.5);  % range 0 to +1, mean = 0.5
%     end
%     
%     if strcmp(stim.carrA.type, 'synthTex')
%         if stim.carrA.nPatt < 1530
%             bubbleMask = carrPattA;
%         end
%     end
%         
%     % adjust carrier mean level, so that AFTER enveloping the resultant mean will then = 0.5:
%     if ~(stim.carrA.mpType==10 || stim.carrA.mpType==11)   % ##-CB  i.e., if NOT a gaussian
%         carrPattA = carrPattA - (sum(sum(carrPattA.*envPatt))-0.5)/sum(sum(envPatt)); 
%         carrPattA = 2*carrPattA;
%         fprintf(1,'adjusting mean of carrPattA\n');  % $$ diag
%     else
%         carrPattA = (carrPattA - 0.5)*2;  % ## -> about same scale (??) but not fancy envelope-correction ?
%         % $$ ## evidently this is has a small but very significant negative offset
%         %           whose magnitude increases with the number of micropatterns
%         %  $$ ## -> was this actually present in original textures, but only noticeable after enveloping ??
%         fprintf(1,'NOT adjusting mean of carrPattA\n');  % $$ diag
%     end
%     
% %     figure(2);
% %     subplot(2,3,4);  imagesc(carrPattA,[-2 2]); title('re-adjusted carrPattA'); colorbar;  % $$ diag
% %     figure(3);
% %     subplot(2,3,4); plot(1:480,carrPattA(240,:), 1:480,0*ones(1,480));    % $$ diag
% 
%     thisStimImage = carrPattA.*envPatt;  % re-envelope, now with zero-mean -> range -.5 to +.5
%     
%     if stim.carrA.whiten
%         thisStimImage = thisStimImage + 0.5;
%         imgWhite = imgWhiten(thisStimImage,stim.carrA.filtSigma,stim.carrA.filtOrder);
%         imgWhite = imgZclip(imgWhite,stim.carrA.sdClip); % clip at e.g. +/-2*SD
%         imgWhite = imgWhite - 0.5;      % range -.5 to +.5
%         imgWhite = winMask.*imgWhite;   % apply cosine-window taper
%         thisStimImage = imgWhite;
%     end
%     
%     thisStimImage = 2*thisStimImage;    % -> range -1 to +1
%     
% %     figure(2);
% %     subplot(2,3,5);  imagesc(thisStimImage,[-1 1]); title('unadjusted thisStimImage'); colorbar;  % $$ diag
% %     % $$ ## - for gabors, imagesc w.o. [-1 1] reveals possible "rogue pixel" ~ -4 ??
% %     figure(3);
% %     subplot(2,3,5); plot(1:480,thisStimImage(240,:), 1:480,0*ones(1,480));   % $$ diag
% 
%     if ~(stim.lumConAdd==0)    % add a specified luminance contrast (note this must come AFTER conScl)
%         lumPatt = stimMakeEnvPattPsy(stim.env,imgSize); % luminance version of envelope, range -1 to +1        
%         thisStimImage = thisStimImage + (stim.lumConAdd/100)*lumPatt.*winMask; 
%     end        
%       
%         
%     thisStimImage = scr.gray + scr.inc*thisStimImage;   % range 0 to 255
%     
% %     figure(2);
% %     subplot(2,3,6);  imagesc(thisStimImage,[0 255]); title('final thisStimImage'); colorbar;  % $$ diag
% %     figure(3);
% %     subplot(2,3,6); plot(1:480,thisStimImage(240,:), 1:480,127.5*ones(1,480));    % $$ diag

    
if strcmp(stim.modType,'quilt') ...
       || strcmp(stim.modType,'quilt_con') % quilts a la Landy-Oruc
    
    envPattA = sqrt((1 + (stim.modDepth/100)*envPatt)/2.0).*winMask;  % range -1 to +1
    envPattB = sqrt((1 - (stim.modDepth/100)*envPatt)/2.0).*winMask;      
    % note that within the disk, (envPattA.^2 + envPattB.^2) always sums to 1.0, seamlessly across the contour
%      
%     if(preComp)
%         carrPattA  = double(carrAPattsGlobal{pattInd});
%         carrPattB  = double(carrBPattsGlobal{pattInd});
%         flipStateA = 'nn';
%         flipStateB = 'nn';
%         rSeedA     = ClockRandSeed;
%         rSeedB     = ClockRandSeed;        
%     else
          
        [carrPattA flipStateA rSeedA] = stimMakeCarrPattPsy(stim.carrA,scr.gray,stim.size); 
                                                                    % sdClipped, range = 0 to zMax (unknown)
        [carrPattB flipStateB rSeedB] = stimMakeCarrPattPsy(stim.carrB,scr.gray,stim.size); 
         
%     end
    
%     if (stim.carrA.mpType==10 || stim.carrA.mpType==11)
%         mpType = 'lum';    % luminance-defined, i.e. pos_ or neg_gaussian
%     else
%         mpType = 'con';     % contrast-defined, e.g. edgelets
%     end
%     
%     if strcmp(mpType,'con')   % ##-CB  i.e., if NOT a gaussian
%         carrPattA = carrPattA - mean(carrPattA(:));     % zero-mean
%     elseif strcmp(mpType,'lum')   % i.e., pos_gaussian or neg_gaussian
%         carrPattA = carrPattA - 127.5;   % (this is correct) 
%     end
    
    carrPattA = carrPattA./max(abs(carrPattA(:))); % range -1 to +1, with mean = 0.0
    carrPattA = (carrPattA + 1.0)/2.;               % range 0 to +1, with mean = 0.5
    
%    fprintf(1,'carrA mean=%f, max=%f, min=%f\n',mean2(carrPattA),max(max(carrPattA)),min(min(carrPattA))); 
% ranges 0 to 1.0  (0-.5 = dark, .5-1.0 = light)

    cv = std2(carrPattA)/mean2(carrPattA);  % "raw" rms-contrast (coefficient of variation), typically 0-1 
    conSclA = stim.carrA.con / (100*cv);    % adjust carrier contrast to specified rms-contrast

    avg = 2*abs(mean2(carrPattA-0.5));      % "raw" mean luminance - max possible = 0.5 ?
    lumSclA = stim.carrA.con / avg;         % adjust carrier mean to this specified value, on scale 0-1
    
%     if strcmp(mpType,'con')   % ##-CB  i.e., if NOT a gaussian
%         carrPattA = 0.5 + conSclA*(carrPattA - 0.5);  % range 0 to +1, mean = 0.5
%         carrPattA = carrPattA - (sum(sum(carrPattA.*envPattA))-0.5)/sum(sum(envPattA));    
%     elseif strcmp(mpType,'lum') 
%         carrPattA = (lumSclA*(carrPattA - 0.5))*2;  % $$ ## approx, not incorporating envelope ... 
%     end
    
%fprintf(1,'adj carrA mean=%f, max=%f, min =%f\n\n',mean2(carrPattA),max(max(carrPattA)),min(min(carrPattA)));  
    
    halfPattA = carrPattA.*envPattA;  % re-envelope, now with zero-mean -> range -.5 to +.5
    halfPattA = 2*halfPattA;
    
%     
%     % now handle carrPattB similarly:
%     if strcmp(mpType,'con')    % ##-CB  i.e., if NOT a gaussian
%         carrPattB = carrPattB - mean(carrPattB(:));     % zero-mean
%     elseif strcmp(mpType,'lum')   % i.e., pos_gaussian or neg_gaussian
%         carrPattB = carrPattB - 127.5;
%     end
    carrPattB = carrPattB./max(abs(carrPattB(:))); % range -1 to +1, with mean = 0.0
    carrPattB = (carrPattB + 1.0)/2.;               % range 0 to +1, with mean = 0.5

%fprintf(1,'carrB mean=%f, max=%f, min=%f\n',mean2(carrPattB),max(max(carrPattB)),min(min(carrPattB))); 
% ranges 0 to 1.0  (0-.5 = dark, .5-1.0 = light)

    cv2 = std2(carrPattB)/mean2(carrPattB);  % i.e., coefficient of variation, typically 0-1        
    conSclB = stim.carrB.con / (100*cv2);
    
    avg = 2*abs(mean2(carrPattB-0.5));     % "raw" mean luminance
    lumSclB = stim.carrB.con / avg; 
      
%     if strcmp(mpType,'con')    % ##-CB  i.e., if NOT a gaussian
%         carrPattB = 0.5 + conSclB*(carrPattB - 0.5);  % range 0 to +1, mean = 0.5
%         carrPattB = carrPattB - (sum(sum(carrPattB.*envPattB))-0.5)/sum(sum(envPattB)); 
%     elseif strcmp(mpType,'lum') 
%         carrPattB = (lumSclB*(carrPattB - 0.5))*2;  % $$ ## approx, not incorporating envelope ... 
%     end
    
%fprintf(1,'adj carrB mean=%f, max=%f, min=%f\n\n',mean2(carrPattB),max(max(carrPattB)),min(min(carrPattB)));  
    
    % finally, join A to B:
    halfPattB = carrPattB.*envPattB;  % re-envelope, now with zero-mean -> range -.5 to +.5
    halfPattB = 2*halfPattB;
    
    thisStimImage = halfPattA + halfPattB;  % quilt, range -1 to +1
       
    
%     
%     if ~(stim.lumConAdd==0)  % add a specified luminance contrast (note this must come AFTER conScl)
%         lumPatt = stimMakeEnvPattPsy(stim.env,imgSize); % luminance version of envelope, range -1 to +1        
%         thisStimImage = thisStimImage + (stim.lumConAdd/100)*lumPatt.*winMask; 
%     end
%     
%     if stim.carrA.whiten
%         thisStimImage = (thisStimImage + 1.0)/2.0;
%         imgWhite = imgWhiten(thisStimImage,stim.carrA.filtSigma,stim.carrA.filtOrder);
%         imgWhite = imgZclip(imgWhite,stim.carrA.sdClip); % clip at e.g. +/-2*SD
%         imgWhite = imgWhite - 0.5;      % range -.5 to +.5
%         imgWhite = winMask.*imgWhite;   % apply cosine-window taper
%         thisStimImage = 2*imgWhite;
%     end
%   
%     if strcmp(stim.modType, 'quilt_con')
%         if strcmp(stim.conLoc, 'A')
%             envPattC = ((stim.conAdd/100)*envPatt + 1.0)/2.0;  % range 0 to 1
%         else
%             envPattC = ((stim.conAdd/100)*envPatt*-1 + 1.0)/2.0;  % range 0 to 1
%         end
%          envPattC = envPattC.*winMask;
%          thisStimImage = envPattC.*thisStimImage;
%     end
    
    
    thisStimImage = scr.gray + scr.inc*thisStimImage;   % range 0 to 255
%     imgPrintStats(thisStimImage);

    envPatt = envPattA; % (for diagnostics)
%     
% elseif strcmp(stim.modType,'lumCarr')   % luminance modulation, but with added uniform carrier
%     error('this needs some work !');  % $$-> impose reforms of 'con'
%     [carrPattA flipStateA rSeedA] = (stim.carrA.con/100)*stimMakeCarrPattPsy(stim.carrA,scr.gray,stim.size);  
%                                                                                               % range -1 to +1
%     envPatt = (stim.modDepth/100)*envPatt;
% 
%     envPatt = envPatt + carrPattA;  % possible range = -2 to +2  ( !! )
% 
%     if (max(max(envPatt))>1.0 || min(min(envPatt))<-1)
%         error('luminance-modulation + carrier too large:  reduce modDepth or carrier contrast !');
%     end
% 
%     envPatt = winMask.*envPatt;     % apply window mask                           
% 
%     thisStimImage = scr.gray + scr.inc *envPatt; % range 0 to 255

else
    error('unrecognized stim.modType !');
end


if diagnostics.verbosity
    fprintf(1,'\nthisStimImage, final version:\n');
    imgPrintStats(thisStimImage);
end        

if (max(thisStimImage(:))>scr.white || min(thisStimImage(:))<scr.black)   % final bounds-check
    thisStimImage = min(thisStimImage,scr.white);
    thisStimImage = max(thisStimImage,scr.black);
%     fprintf(1,'\nafter clipping to enforce pixel-depth limitations:\n');
%     imgPrintStats(thisStimImage);
end

    
%  diagnostic plot of image components
if diagnostics.images    
    [scr.console.xSize, scr.console.ySize] = Screen('WindowSize',0);  % size of console screen    
    yOff = 50;  ySize = round(.9*scr.console.ySize);   % (we want it really tall)
    xOff = 10;  xSize = round(.5*ySize);    % proportionate for 4x2    
    fig.pos    = [xOff yOff xSize ySize];
    fig.handle.imgParts = figure('position',fig.pos,'name','image diagnostics');
    figure(fig.handle.imgParts);
    
    if exist('carrPattA','var')
        subplot(4,2,1); imagesc(carrPattA); axis image; colorbar; colormap(gray); title('carrPattA');

        if strcmp(stim.modType,'quilt')
            subplot(4,2,2); imagesc(carrPattB); axis image; colorbar; title('carrPattB');
        end
    end
    subplot(4,2,3); imagesc(envPatt); axis image; colormap(gray); colorbar;  title('envPatt');
    
    subplot(4,2,4); imagesc(thisStimImage,[0 255]); axis image;  colorbar;  title('thisStimImage');

    filtSize = round(stim.size/5)+1;  % must be ODD, for entropyMap or stdMap
    
    filtLow = fspecial('gaussian',filtSize,filtSize/3);  
    
    lumMap = filter2(filtLow,thisStimImage-scr.gray,'same');
    subplot(4,2,5);  imagesc(lumMap); axis image;  colorbar;  title('lumMap');
    
%    rmsMap = stdfilt(thisStimImage,ones(filtSize,filtSize));
    rmsMap = (thisStimImage-mean2(thisStimImage)).^2;  % square each value
    rmsMap = filter2(filtLow,rmsMap,'same');   % locally average these
    rmsMap = sqrt(rmsMap);
    subplot(4,2,6);  imagesc(rmsMap); axis image; colorbar;  title('rmsMap');
      
    drawnow;
end



% if (diagnostics.halves==1)  % NOTE:  somehow (???) this has gotten zeroed
% before getting here !
    extMat=zeros(stim.size,stim.size);
    for i=1:stim.size
        for j=1:stim.size
            if j>=i 
                extMat(i,j)=1;
            else
                extMat(i,j)=0;
            end
        end
    end
    
    %$$ Warning: Obsolete syntax. In future versions IMROTATE will return the result in ans 
    %       instead of displaying it in figure.
    if stim.env.ori==45
        %imrotate(extMat,90);
    end
    
    half1=extMat.*thisStimImage;
    [r c half1]=find(half1);
    extMat=~extMat;
    half2=extMat.*thisStimImage;
    [r c half2]=find(half2);
    meanhalf1= mean2(half1);
    meanhalf2 = mean2(half2);
    rms1 = sqrt(mean2((half1-meanhalf1).^2));     % rms contrast
    rms1 = rms1 / meanhalf1;     % this should now be dimensionless
    rmsCon1 = rms1 * 100;  % express as a percentage
    %disp(rmsCon1)
    rms2 = sqrt(mean2((half2-meanhalf2).^2));     % rms contrast
    rms2 = rms2 / meanhalf2;     % this should now be dimensionless
    rmsCon2 = rms2 * 100;  % express as a percentage
    %disp(rmsCon2)
    %imgPrintStats(thisStimImage);
    

%     fprintf(1,'%5.5f %5.2f %5.2f %5.2f %5.2f\n',stim.modDepth, meanhalf1, meanhalf2, rmsCon1, rmsCon2);
  
    lumMcon = 100*(meanhalf1 - meanhalf2) / (meanhalf1 + meanhalf2);
    rmsMcon = 100*(rmsCon1   - rmsCon2)   / (rmsCon1   + rmsCon2);
    
 %   fprintf(1,'%5.5f %5.2f %5.2f %5.2f    %5.2f %5.2f %5.2f\n',stim.modDepth, meanhalf1, meanhalf2, lumMcon,   rmsCon1, rmsCon2, rmsMcon);
    
 
    
    %fprintf(1,'modDepth = %5.1f\n', stim.modDepth);
    %fprintf(1,'Mean of first half = %5.1f , RMS Contrast of First half is %5.1f\n',meanhalf1,rmsCon1);
    %fprintf(1,'Mean of Second half = %5.1f , RMS Contrast of Second half is %5.1f\n',meanhalf2,rmsCon2);
% end

