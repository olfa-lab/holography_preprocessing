function bg_sub_mov = morph_bg_sub(mov, morph_radius, blur_sigma)
% Use blurring + morphological opening to get the background of a movie and
% subtract it
%   morph_radius: morphological disk element radius (choose such that all
%   elements of the foreground are within this distance from a background
%   space)
%   blur_sigma : sigma for initial gaussian spatial blurring to reduce
%   noise prior to morphological operations

nframes = size(mov,3);

se = strel('disk', morph_radius); % filtering element

for iframe=1:nframes
   initframe = mov(:,:,iframe);
   frame = imgaussfilt(initframe, blur_sigma);
   frame = imopen(frame,se);
   mov(:,:,iframe) = initframe - frame;
end

bg_sub_mov = mov;

end

