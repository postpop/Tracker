function fp = fixOrientations(fp)
% fp = fixOrientations(fp)
%
% ARGUMENTS:
%     fp  - FlyPursuit object 
%           or structure with fields TRACKS, NFLIES, ORIENTATION

% 0. extract params
posX = mapFun(@smooth, fp.tracks(:,:,1),10);
posY = mapFun(@smooth, fp.tracks(:,:,2),10);
velX = [zeros(1, fp.nFlies); diff(posX)];
velY = [zeros(1, fp.nFlies); diff(posY)];
angVel = cart2pol(velX, velY);
spd = sqrt(velX.^2 + velY.^2);

oriX = fp.orientation(:,:,1);
oriY = fp.orientation(:,:,2);
angOri = cart2pol(oriX, oriY);
dAngOri = [zeros(1, fp.nFlies); diff(angOri)];

% unwrap angles
dAngOri(dAngOri<-pi) = dAngOri(dAngOri<-pi) + 2*pi;
dAngOri(dAngOri>pi) = dAngOri(dAngOri>pi) - 2*pi;

%% 1. correct random flips (>85degrees rotation within a frame)
dAngOri(dAngOri<-.9*pi) = dAngOri(dAngOri<-.9*pi) + pi;
dAngOri(dAngOri>.9*pi) = dAngOri(dAngOri>.9*pi) - pi;

[oriX, oriY] = pol2cart(cumsum(dAngOri),1);

%% 2. fix simple inversion:
% assume that flies move forwards most of the time - hence orientation and
% velocity vector should be aligned most of the time
% calc velocity and flip if neg.
% prj onto orientation vector to get forward component
for fly = 1:fp.nFlies
   ori = [oriX(:,fly) oriX(:,fly)];
   vel = [velX(:,fly) velX(:,fly)];
   tmp = ori'*vel;% avg. 'angle' between orientation and speed vector
   totalDir(fly) = tmp(1);
end
oriX = bsxfun(@times, oriX, sign(totalDir));
oriY = bsxfun(@times, oriY, sign(totalDir));

%%  
fp.orientation(:,:,1) = oriX;
fp.orientation(:,:,2) = oriY;