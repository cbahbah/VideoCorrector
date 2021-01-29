# VideoCorrector
Create initial video from corrupted video sequence 
by Chahrazade Bahbah (https://github.com/cbahbah/)

## External Libraries
 - OpenCV : https://en.wikipedia.org/wiki/OpenCV


## Methodology
 - First, we extract the frames from the video using openCV methods. 
 - Then, we delete the "outlier" frames. To do so, two methods were tested, **PSNR** and **SSIM**. After several tests, SSIM have been selected since it gives better results. However, defining a treshold to determine if a frame should be kept or not is not straightforward. Several thresholds were tested : average value, median value or standard deviation. This step gives us a list of the unordered "good" frames. The resulting video is **filtered.mp4**. The computation of the SSIM index is heavy in terms of computational time, thus the values are stored in text files (SSIM.txt, meanSSIM.txt) and processed later.
 - To reoder the selected frames, we recompute the SSIM index (SSIMfiltered.txt, meanSSIMfiltered.txt). This time, instead of comparing the mean values, we compare the values frame by frame. For now, we consider that the fist frame of the video is the correct first one. So we start from the first, and move forward. We compare the selected frame with rest to find the closest one. The process is reiterated till all frames have been processed. We create the final video from the ordered frames (**finalVideo.mp4**).


## Future Improvements 
 - The algorithm SURF should be tested to reaorder the frames. It consists in detecting interesting points in each frame and compute the correspondances with the remaining frames. Unlike the SSIM method that compares the luminance and the contrast between the images and output a metric, the SURF method is a local descriptor that give us a number of points correspondances between the frames. Note that this could correct the video sequence if unwanted frames are still remaining. The method was implemented but not tested yet (no time :( ).
 - Also, a method that compute the first or last frame of the video sequence should be implemented. A possible algorithm could be to select the first frame as the one who has the least point correspondances with others. In fact, if a frame is located in the middle of the video sequence, the higher the correspondance is...



## Sources
 - PSNR : https://fr.wikipedia.org/wiki/Peak_Signal_to_Noise_Ratio
 - SSIM : https://fr.wikipedia.org/wiki/Structural_Similarity
 - SURF : https://fr.wikipedia.org/wiki/Speeded_Up_Robust_Features




