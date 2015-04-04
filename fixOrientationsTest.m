cc()
% trunk = '/Volumes/murthy/jan/courtship/dat/141202_1137/141202_1137';%'dat/140827_1329/140827_1329';
trunk = '/Volumes/murthy/jan/courtship/dat/140728_1732/140728_1732';%'dat/140827_1329/140827_1329';
load([trunk '_res']);
fp.vr = VideoReader([trunk '.mp4']);
fp.tracks = p.tracks;
fp.orientation = p.orientation;
%% fix orientations
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

% correct flips (>85degrees rotation within a frame)
dAngOri(dAngOri<-.9*pi) = dAngOri(dAngOri<-.9*pi) + pi;
dAngOri(dAngOri>.9*pi) = dAngOri(dAngOri>.9*pi) - pi;

[oriX, oriY] = pol2cart(cumsum(dAngOri),1);

dAngOriAngVel = angOri - angVel;

% 1. align with speed vector - NO
% idx = abs(dAngOriAngVel)>pi/2 & abs(dAngOriAngVel)<3*pi/2;
% idx = idx & spd>1;
% oriX(idx) = -oriX(idx);
% oriY(idx) = -oriY(idx);

% %% 2. fix flips due to degeneracy of ellipse orientation info:
% clf
% subplot(211)
% plot(dAngOri(2000:3000,1),'.-')
% hold on
% for fly = 1:fp.nFlies
%    jumpIdx = find(abs(dAngOri(:,fly))>.8*pi);
%    jumpIdx(jumpIdx>3000) = [];
%    for i = 1:length(jumpIdx)-1
%       oriX(jumpIdx(i):end,fly) = -oriX(jumpIdx(i):end,fly);
%       oriY(jumpIdx(i):end,fly) = -oriY(jumpIdx(i):end,fly);
%    end
% end
% 
% angOri = cart2pol(oriX, oriY);
% [oriX, oriY] = pol2cart(angOri,1);
% dAngOri = [zeros(1, fp.nFlies); diff(angOri)];
% 
% dAngOri(dAngOri<-pi) = dAngOri(dAngOri<-pi) + 2*pi; % unwrap??
% dAngOri(dAngOri>pi) = dAngOri(dAngOri>pi) - 2*pi; % unwrap??
% 
% subplot(212)
% plot(dAngOri(2000:3000,1),'.-')
% axis(gcas,'tight')
%% 2. fix simple inversion:
% assume that flies move forwards most of the time - hence orientation and
% velocity vector should be aligned most of the time
% calc velocity and flip if neg.
% prj onto orientation vector to get forward component
for fly = 1%:fp.nFlies
   ori = [oriX(:,fly) oriX(:,fly)];
   vel = [velX(:,fly) velX(:,fly)];
   tmp = ori'*vel;
   totalDir(fly) = tmp(1);
end
totalDir
oriX = bsxfun(@times, oriX, sign(totalDir));
oriY = bsxfun(@times, oriY, sign(totalDir));

fp.orientation(:,:,1) = oriX;
fp.orientation(:,:,2) = oriY;


%%

clf
colormap(gray)
fp.initFrame = 20000;
fp.currentFrameIdx = fp.initFrame-1;
history = 12;
cmap = limit(jet(fp.nFlies+4)-.2);
lag = 8;12
for fr = fp.initFrame:fp.NumberOfFrames
    fp.getNextFrame();
    imagesc(fp.currentFrame);
   hold on;
   t = fr+lag;
   plot(fp.tracks(t-(1:history),:,1), fp.tracks(t+-(1:history),:,2),'LineWidth',2);
   flyX = fp.tracks(t,:,1);
   flyY = fp.tracks(t,:,2);
   plot(flyX, flyY,'.','MarkerSize',18);
   plot([flyX;flyX + 50*fp.orientation(t,:,1)], [flyY;flyY + 50*fp.orientation(t,:,2)])
   %plot([flyX;flyX + 50*oriX(t,:)], [flyY;flyY + 50*oriY(t,:)])
   hold off
   title(t)
   set(gca,'XLim',[0 fp.w], 'YLim',[0 fp.h])
   set(gca,'DataAspectRatio',[1 1 1])
   drawnow
end