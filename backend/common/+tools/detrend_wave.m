function wave = detrend_wave( wave, filt, roi )
% Removes the trendline for an input wave
% version 0.1,  20 Feb 2020  Calvin Eiber

if length(filt) == size(wave,1)
  
  time = filt; 
  k = flattopwin( round(20/mean(diff(time))) + 1 ); % 500 Hz ? 
  k = k/sum(k);
  
else k = filt;   
end

for ty = 1:numel( wave(1,1,:) ) % may be 3D or 4D

  U = mean(wave(:,:,ty),2); % mean across electrodes

  U_lf = conv(U,k,'same');
  U_ok = conv(conv(0*U+1,k,'same'),k,'same').^2;
  U_lf = U.*(1-U_ok) + U_lf.*U_ok; 
  
  if ~any(U_lf(roi)), continue, end % zero wave

  le = U_lf(roi) \ wave(roi,:,ty); % Linear scaling term
  wn = wave(:,:,ty) - U_lf(:) * le;
  w0 = conv2(wn,k,'same'); 

  wave(:,:,ty) = wn - w0;

end