function frames = extractFlyFrames(p)
% return frames cropped around the fly and rotated to align fly

p.framesToRead = unique(p.framesToRead);
frames = zeros(length(p.boxW), length(p.boxH), p.imageChannels, p.NumFramesToRead, 'uint8');

flyIdxIdx = p.framesToRead - p.fp.initFrame + 1;% convert frameidx to array indices by correcting for initFrame
if p.flyBoxVidExists
   framesToRead = p.framesToRead - p.fp.initFrame + 1;
else
   framesToRead = p.framesToRead;
end
for f = 1:length(p.framesToRead)
   frame = p.fp.vr.read(framesToRead(f));
   if p.flyBoxVidExists
      % leave frame as is
      flyFrame = frame;
   else
      % convert track idx (in cropped frame) to global frame idx
      flyX = round(p.fp.tracks(framesToRead(f),:,1) + min(p.fp.boundsX));
      flyY = round(p.fp.tracks(framesToRead(f),:,2) + min(p.fp.boundsY));
      flyFrame = fix.extractFlyBox(frame,flyX,flyY,p.boxW,p.boxH);
   end
   singleFlyFrame = flyFrame(:, (p.flyIdx(flyIdxIdx(f),p.currentFly)-1)*length(p.boxH) + (1:length(p.boxH)) ,:);
   
   % rotate frames according to fly orientation
   if p.flyBoxVidExists
      thisAng = p.flyAngle(framesToRead(f)+p.fp.initFrame -1, p.flyIdx(flyIdxIdx(f),p.currentFly));
   else
      thisAng = p.flyAngle(framesToRead(f), p.flyIdx(flyIdxIdx(f),p.currentFly));
   end
   singleFlyFrame = imrotate(singleFlyFrame, thisAng+90, 'crop');
   
   frames(:,:,:,f) = singleFlyFrame;
end
