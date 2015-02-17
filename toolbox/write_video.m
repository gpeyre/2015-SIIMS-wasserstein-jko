function write_video(A, filename, filetype, options)

% write_video(A, filename, filetype, options);

switch filetype
    case 'gif'
        q = 128; % color map size
        if size(A,4)>1
            % conversion to indexed
            X = [];
            [X(:,:,1,1),map] = rgb2ind(A(:,:,:,1),q);
            for i=2:size(A,4)
                [X(:,:,1,i),map] = rgb2ind(A(:,:,:,i),map);
            end
            imwrite(X+1, map, [filename '.gif'], 'DelayTime',1/20, 'loopcount', Inf);
        else
            res = @(A)reshape(A, [size(A,1) size(A,2) 1 size(A,3)]);
            CM = gray(q);
            imwrite(rescale(res(A),1,q), CM, [filename '.gif'], 'DelayTime',1/20, 'loopcount', Inf);
        end
    case {'mp4' 'avi'}
        vidObj = VideoWriter( filename);
        vidObj.Quality = 30;
        vidObj.Quality = 50;
        % vidObj.VideoFormat = 'Grayscale';
        %vidObj.CompressionRatio = 20;
        open(vidObj);
        if size(A,4)>1
            for i=1:size(A,4)
                writeVideo(vidObj,A(:,:,:,i));
            end            
        else
            for i=1:size(A,3)
                writeVideo(vidObj,A(:,:,i));
            end
        end
        close(vidObj);
end