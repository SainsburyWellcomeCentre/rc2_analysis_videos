% visualize_face_and_body_movement.m
%
% Description:
%   This script selects specific frames from video recordings of face and body 
%   movements to visualize motion. It extracts two consecutive frames from each 
%   specified point in the video to highlight changes in pixel intensity, 
%   allowing for a clear view of movement areas in face and body recordings.
%
% Workflow:
%   1. **Video Loading and Frame Selection**:
%      - Loads two video files representing face (`camera0`) and body (`camera1`) 
%        movement recordings.
%      - For each video, selects a specific frame based on `selected_frame_number`.
%      - Displays the selected frame with annotations of the current time and frame 
%        number.
%
%   2. **Motion Visualization**:
%      - Captures two consecutive frames and calculates the pixel-wise difference, 
%        emphasizing regions of movement.
%      - Displays the difference between frames as a heatmap to visualize movement 
%        intensity.
%      - Uses a custom color map to highlight movement areas.
%
% Usage:
%   Specify `selected_frame_number` for the frames to visualize, and update 
%   `path_to_example_videos` with the paths to the video files. Run the script to 
%   generate visualizations showing regions of face and body movement.


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