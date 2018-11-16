# ssp7
## SSP Course Project Group 7
Aarsh Prakash Agarwal 150004 aarshp@iitk.ac.in
Harshvardhan 150283 harshv@iitk.ac.in
Shubhanshu Khandelwal 150705 shubhke@iitk.ac.in
Krishna Chandra 150337 krishc@iitk.ac.in
Vishal Srivastava 150820 vishalsr@iitk.ac.in


## Files for Paper-2
### sampler.m: 
This is the main script for the second paper's implementation. It defines the 3 parts of the algorithm, source motion, sampling possible source positions and updating the weights for the sampled source position using SSP\_EM function on the convolution of the RIR and STFT of the observed speech signal. Then it plots the KDE contour plots for each time step.

### stft.m: 
Function takes in the time domain signal values, sampling rate, window, and hop length(here it is the same as window length) and outputs the STFT of the given signal. This function has been ported from MATLAB File Exchange.

### RIR: 
Function has been ported from the RIR library. It computes the RIR for given source-microphone geometries when the  reverberation time and room dimensions have been supplied. 

### SSP_EM.m:
This file runs the subroutine EM algorithm for at max 50 iterations. It uses certain threshold for the convergence. This file assumes following input from file sampler.m
    1. microphone position
    2. sampled state(x,y,v_x,v_y) of particles
    3. STFT of the signal recieved by microphones.(size = 2 X K; K = no. of frequency bins.)
    4. $\epsilon$ to calculate the Gamma(t,k) for each instant t and each frequency bin k
This subroutine returns probability vector which inturn is used to update weights w_j^(t). This subroutine calls complex_gauss and calc_distance as helper functions.

### Complex_gauss.m
Function that takes point and covariance matrix (as a complex number) as input and returns the probability of the point from the circular Gaussian distribution with zero mean and given covariance matrix.

### calc_distance.m
Function that takes two matrices as input and returns the maximum of the absolute value of the difference between the two.

## Files for paper-1
There are basically two files. The first file implements the algorithm and plots the spectrum while comparing its with state of the art MUSIC and MVDR algorithms, while the second file implements the Monte Carlo simulation and compares their performance on different SNR values. Following is the description:
### doa.m
The file implements the signal model and IG algorithm, and compares its spatial spectrum with MUSIC and MVDR. MUSIC and MVDR are implemented via in-built phased library functions present in MATLAB.
### monte_carlo_doa.m
The file implements the Monte Carlo simulation for all the three algorithms. The simulation is run for different values of SNR depending upon the problem. Sources are kept to 2 in this file.

