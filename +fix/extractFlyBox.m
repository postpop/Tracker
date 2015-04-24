function flyFrame = extractFlyBox(frame, flyX, flyY, flyBoxWidthIdx, flyBoxHeightIdx)
% convert track idx (in cropped frame) to global frame idx
[width, height, colorChannels] = size(frame);
nFlies = size(flyX,2);
flyFrame = zeros(length(flyBoxWidthIdx), nFlies*length(flyBoxHeightIdx), colorChannels,'uint8');
for fly = 1:nFlies
    flyBoxX = limit(flyX(fly) + flyBoxWidthIdx, 1, width);
    flyBoxY = limit(flyY(fly) + flyBoxHeightIdx, 1, height);
    flyFrame(1:length(flyBoxX), (fly-1)*length(flyBoxWidthIdx) + (1:length(flyBoxY)),:) = frame(flyBoxX, flyBoxY,:);
end

