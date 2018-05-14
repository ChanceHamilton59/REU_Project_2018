% FUNCTION: drawSynthTexNatCD - draw a single micropattern texture
%  
% usage:  littletex = drawSynthTexNat(fTSize, npatt, mpCode, mpOris, mpSizes)
%
% args:  fTsize -- final texture size (480)
%        nPatt  -- number of micropatterns (some multiple of 85)
%        mpCode -- numerical index for micropattern type (same as mpCode)
%                       1 - 2 frequencies, intact
%                       2 - broadband,     intact ("edgelet")
%                       3 - gabor,         intact
%                       4 - 2 frequencies, phase randomized
%                       5 - broadband,     phase randomized
%                       6 - gabor,         phase randomized
%                       7 - 2 frequencies, orientation randomized
%                       8 - broadband,     orientation randomized
%                       9 - gabor,         orientation randomized (trivial)
%                       10- pos_gaussian   intact ("light_blob")
%                       11- neg_gaussian   intact ("dark_blob", negative of #10)
%        mpOris  -- some subset of [1:30:360]
%        mpSizes -- 4 element vector specifying sizes of micro-patterns
%                       In Liz's experiments this was [16 32 64 128]
%                       To present one size only pass in a vector with all
%                       entries equal (e.g. [32 32 32 32]).
%
%
% global:
%       mpLibrary - micropattern library (produced by makeEdgeletLib_os)
%       s - index for which size of micropattern
%       o - index for which orientation of micropattern
%       nVari - number of variations of scrambled micropatterns in library 
%
% textures are [0-255] with mean/background = 127.5
% textures are NOT rms-normalized; and for gaussian textures, not luminance-equated

% written by Elizabeth Arsenault (now Zavitz), 2012-2013


% example usage:   
% global mpLibrary s o nVari;
% load mpLibrary_os;
% fTSize = 480;
% npatt = 85*7;
% mpCode = 9;  % gabors
% mpOris = o;  % mix of all orientations
% tex = drawSynthTexNat(fTSize, npatt, mpCode, mpOris);
% imagesc(tex);  colormap(gray); axis image;

% TO DO:
% in conditional near end, avoiding removal of mean, avoid dependence on mpCode
%       to make this robust for different libraries

% recent changes:
% 16 Sept 2013, CB: changed RandSample to randsample
% 20 Sept 2013, CB: modify final "normalization" to ensure zero-background for gaussians
% 
% CD is Chris DiMattina, visitor from FGCU (cdimattina@fgcu.edu)
% 23 May  2016, CD: added ability to handle single micro-pattern sizes 
%
%

% needs to run:
% (PsychToolbox - for RandSample, which is NOT same as Matlab's randsample !)

% see also:
%   test_drawSynthTexNat.m 


function littletex = TextureGenCH(fTSize, npatt, mpCode, mpOris, mpSizes)

% mpLibrary variables:
% mpLibrary -- the library itself
% s         -- index for mp by size
% o         -- index for mp by orientation
% nVari     -- number of variations of scrambled micropatterns in library
global mpLibrary s o nVari;

%load mpLibrary_os

tSize = fTSize + 32;              % texture size

% CD: 05-24-2016
% Commented this out - now passed in as parameter
% mpSizes = [16 32 64 128];       % micropattern sizes

texture = zeros(tSize); %initialize texture canvas

% build array of sizes that will be 1/f

% rs = [repmat(mpSizes(1),1,round((64/85)*npatt)) ...
%       repmat(mpSizes(2),1,round((16/85)*npatt)) ... 
%       repmat(mpSizes(3),1,round((4/85)*npatt)) ... 
%       repmat(mpSizes(4),1,round((1/85)*npatt))];

% rs = Shuffle(rs);  %so that the sizes are thrown down in random order

rs = PsychRandSample(mpSizes, [1 npatt]);

ro = PsychRandSample(mpOris, [1 npatt]); %random orientations for micropatterns
% note this is PsychToolbox function, NOT same as Matlab's randsample.m 

whichV = 1; % which version of the micropattern to use. 
            % Randomized only if phase or orientation scrambled
            
% this loop places one micropattern on the canvas with each iteration    

% Degugger area disp(s)

for i = 1:npatt 
    
    
    if mpCode > 3 && mpCode <= 9 %if necessary, pick a random mp version
        whichV = round(1 + (nVari-1)*rand(1));
    end

    
    si = find(s == rs(i)); %obtain library indicies for size and orientation
    oi = o == ro(i);

  
    
    
    micropattern = mpLibrary{mpCode, si, oi, whichV}; %look up micropattern
    
    posx = round(1 + (length(texture)-rs(i))*rand(1));%randomize position
    posy = round(1 + (length(texture)-rs(i))*rand(1));
    
    texture(posx:posx+rs(i)-1, posy:posy+rs(i)-1) = ...  %place micropattern in position   
        texture(posx:posx+rs(i)-1, posy:posy+rs(i)-1)+micropattern;   
end

% trim and window
trim = (tSize-fTSize)/2;
littletex = texture(trim:end-trim-1, trim:end-trim-1);

% normalize to 0-255, with mean/background=127.5
if mpCode<10   % i.e. not a gaussian, so everything should be zero-balanced
     littletex = littletex - mean2(littletex);         % remove mean
end

littletex = littletex./max(max(abs(littletex)));  % range between -1 and +1

littletex = 127.5 + 127.5*littletex;                % range 0-255
% i.e. for pos_gaussian or neg_gaussian, ensure "background"=127.5


 
 