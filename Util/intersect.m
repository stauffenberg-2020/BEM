function vq = intersect(x,v,xq)
    % Note : This script does not extrapolate the curves to find the
    % intersection point, Update may be required
    
    y1 = x;
    x = v;
    y2 = repmat(xq,length(x),1);

    dy1y2 = y2-y1;                                                  % Subtract To Create Zero-Crossings
    idx = find(diff(sign(dy1y2)));                                  % Find Approximate Indices Of Zero-Crossings
    idx = idx(1); % Making sure it is a scalar for coder compatibility
    vq = zeros(1,idx); % or xv
    yv = zeros(1,idx);
    for k = 1:numel(idx)
        idxrng = max(idx(k)-2,1) : min(numel(x),idx(k)+2);          % Restrict Index Range To Meet Requirements For ‘interp1’
        vq(k) = interp1(dy1y2(idxrng), x(idxrng), 0);               % Calculate ‘x’ Coordinates Of Zero-Crossings
        yv(k) = interp1(x(idxrng), y2(idxrng), vq(k));              % Calculate ‘x’ Coordinates Of Zero-Crossings (Should Be Close To Zero For Most, Although Not All)
    end
    
%     plot(y1, x,'--k*');hold on;plot(y2, x,'--b.');plot(yv,vq,'r*');
%     xlabel('x'); ylabel('v');
    
end