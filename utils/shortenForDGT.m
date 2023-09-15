function [firstIdx, L] = shortenForDGT(w, a, s, f, offset)

    firstWin = ceil((s - ceil(w/2))/a) + 1;
    firstWinPos = 1 + (firstWin - 1)*a;
    firstWinPos = firstWinPos + offset;

    if firstWinPos - a + ceil(w/2) >= s + 1
        firstWin = firstWin - 1;
        firstWinPos = firstWinPos - a;
    end
    
    firstIdx = firstWinPos - ceil(floor(w/2)/a)*a;

    lastWin = firstWin + floor((f + floor(w/2) - firstWinPos)/a);
    lastWinPos = firstWinPos + (lastWin - firstWin)*a;

    lastIdx = lastWinPos + ceil(floor(w/2)/a)*a;
    
    L = lastIdx - firstIdx + 1;

end