# ECG-R-peaks-detection
R peaks detection in ECGs using wavelet decomposition and higher statistics, implemented in MATLAB

Apostolou Orestis, 26/04/2020

This project is implemented in MATLAB. The wfdb toolbox must be installed in order to execute this. You can find instructions here: https://archive.physionet.org/physiotools/matlab/wfdb-app-matlab/

The task of the project is to identify the R peaks of an ECG.
People can easily detect R peaks in a cardiogram, as they are usually some local maxima of the ECGs. The task here is to create a completely automated algorithm to do this job.

This task may be easy for normal cardiograms, but can prove really challenging for patients with arrythmia.

All the data used, were collected by MIT Arrhythmia Database (https://www.physionet.org/content/mitdb/1.0.0/).
In the project, only 5 cases were used, but it's up to you to modify the code to use it anywhere you need.

This is how the algorithm works:

First of all, it reads the doctors' annotations - the ground truth. Then it reads the complete (digital) ECG. After that, it applies a discrete wavelet decomposition to the waveform. In this project, we chose to apply MATLAB's fk4 wavelet decomposition, and to use the 2nd level's details. By using a level of details in our analysis, we want to get rid of the lower frequencies that exist. We don't want to use the first level of details, as it consists mainly of recording noise.

After aqcuiring the 2nd level of details, our task is to collect all the local maxima that reside above a threshold that we must find. Again, a supervising person can easily find the appropriate threshold by looking the waveform, but we need to have a completely automated algorithm.

Thus, we used a statistic metrics called percentile. (https://en.wikipedia.org/wiki/Percentile)
Specifically, we used the 98th percentile of the full waveform as our threshold. (98th was used after tests)

The results of this technique were quite satisfying. Using a 4 maximum samples margin between the ground and predicted truth (sampling frequency was 360 samples per second), the precision metrics for the 5 patients that we examined, were between 90% -99%; the recall metrics were between 89%-99%. 



Files in the repo:

ECGanalysis.m: Run this script without any arguments, in order to make R peak detection for case 106. You can see many plots and evaluation metrics. Use this script to understand how the algorithm works.

fullECGanalysis.m: Run this script without any arguments, to view the evaluation metrics for all 5 cases examined.

evaluation.m: The script that implements evaluation.

thresholdAdjustment.m: The script that indicates how the tests led us to use 98th percentile.

ECGdata.zip: The 5 cases we examined compressed. You must uncompress it in the same folder as the scripts above.



Feel free to clone, download or use this project for your research, or contribute for something that you would like to add.

For any questions contact me at orestisapostolou@yahoo.gr

