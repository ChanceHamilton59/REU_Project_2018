function envPatt = stimMakeEnvPattPsy(env,imgSize)
%
%stimMakeEnvPattPsy  - return carrier pattern, for use in psychophysics
%  usage:  envPatt = stimMakeEnvPattPsy(thisStim.env,imgSize)
%
%	env  = structure containing stimulus parameters:
%       env.type = 'gabor', 'halfDisc', 'gaussian'
%       env.taper = cosine-taper fraction (0.0 to 1.0), for 'halfDisc'
%       env.lambda, .aspect, .spBand, ... - parameters for 'gabor'
%       env.phase MUST be 0 or 180
%       env.shift = shift of boundary from centre, in pixels, for 'halfDisc'
%   imgSize = image matrix size (must be square), in pixels
%
%   returns envPatt, with values ranging from -1 to +1

% see also:
%   stimMakeCarrPattPsy.m

% RECENT CHANGES:
%   25 Oct 07, CLB: for 'halfDisc', invert thisOri for consistency with response keys,'gabor'
%   26 Oct 07, AY:  'halfDisc' with cosine-taper
%   29 Oct 07, CB:  eliminate thisOri, thisPhase;  replace size by imgSize
%   2  Dec 07, CB:  env.shift - shifted boundary, for halfDisc

if strcmp(env.type,'gabor')       
    a = cos(env.ori*pi/180) * 2*pi/env.lambda;
    b = sin(env.ori*pi/180) * 2*pi/env.lambda;      
    [x y] = meshgrid(-imgSize/2:imgSize/2-1,-imgSize/2:imgSize/2-1); 
                                        % 256^2 matrices, values range -128 to +128
    gratingMatrix = sin(a*x + b*y + env.phase*pi/180);                     
    envSigCross = bw2sigma(env.spBand,env.lambda); % sigma along orthogonal axis
    envSigAxial = envSigCross * env.aspect;        % sigma along axis of orientation
    sigX = envSigCross;     sigY = envSigAxial;     
    gaussMask = exp(-((x.^2)/(2*sigX^2)) -((y.^2)/(2*sigY^2))); % 2-d gaussian, range 0-1       
    envPatt = gratingMatrix .* gaussMask;  % Gabor:  values range -1 to +1
    
elseif strcmp(env.type,'halfDisc')
    tw=env.taper*imgSize;        % define taper width
    if isfield(env,'shift');
        bf = env.shift;     % shift of boundary from centre, in pixels
    else
        bf = 0;
    end
        
    p=zeros(1,imgSize);
    for i=1:imgSize
         if i<= round(bf+imgSize/2-tw/2)
            p(i)=1;
         elseif i>= round(bf+imgSize/2+tw/2)
            p(i)= -1;
        else 
             p(i)=cos((i-round(bf+imgSize/2-tw/2)) *pi/tw);  % create 1-D cosine taper
        end
    end
    
    envPatt = repmat(p,imgSize,1);     % make it 2d
            
    if env.phase==180
        envPatt = -1.*envPatt;    % if phase is 180 inverse
    elseif env.phase==0
        
    else
        error('invalid value of "env.phase"');
    end
        
    envPatt = imrotate(envPatt,-env.ori,'nearest','crop');    % rotate to env.ori
    
%     fprintf(1,'\n\nenv.phase = %f, env.ori = %f\n\n',env.phase,env.ori);  % diag
%     
%     figure(10);  imagesc(envPatt);  axis image;  axis off; colormap(gray); colormap; % diag
%     title(sprintf('env.phase = %f, env.ori = %f\n\n',env.phase,env.ori));
%     drawnow
    
else
    error('unrecognized env.type');
end

end