% Choose the frame you would like to use to visualize face and body
% movement

selected_frame_number = [60 1230];
for i = 0:1
    path_to_example_videos = sprintf('E:\\mvelez\\mateoData_cameras\\CA_176_1_rec1\\camera%d.mp4', i);

    video =  VideoReader(path_to_example_videos);
    

    figure(i + 1);
    subplot(1, 2, 1);
    counter = 0;
    while counter < selected_frame_number(i + 1)
        frame = readFrame(video);
        imshow(frame)
        title(sprintf("Current Time = %.3f sec, frame: %.0f",video.CurrentTime, counter))
        pause(2/video.FrameRate)
        counter = counter + 1;
    end
    
    
    subplot(1, 2, 2);
    frame_first = frame(:, :, :);
    frame = readFrame(video);
    frame_second = frame(:, :, :);

    hold on
    delta = double(frame_first - frame_second);
%     delta(delta(:, :, 1) < 3) = nan;
%     hb = imshow(delta(:, :, 1) * 20, Colormap = copper,Interpolation="bilinear")
    custom_map = [1 1 1
                  1 0 0];
    hb = imshow(delta(:,:, 1), Colormap = custom_map)
%     hb.AlphaData = 0.2;
end