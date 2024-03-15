FROM EMAIL DATED 12/09/2023 from Marieke Mur (WU)

Hi Andrew,

I have, thanks for the nudge! I'm attaching MATLAB code for computing and visualizing RDMs for the monkey ephys data. 

I suggest running START_A, START_B, START_C4a, START_D4, START_D5. This will give you RDM and MDS movies, using the spike rate distance (SRD) as a dissimilarity measure, which we developed for this project some years ago.

I've also included prior versions of START_C and START_D, which compute different dissimilarity measures, but we landed on the SRD because it has a straightforward interpretation (difference in spike rate between two stimuli, averaged across neurons), it is cross validated (meaningful zero point), and does not require multivariate noise normalization (which is great but requires simultaneous measurement of neurons to estimate their covariance matrix). Those prior versions of the code may need some debugging.

I'm also attaching the RSA toolbox code called from the main START scripts. Most START scripts should be documented well enough to figure out what the code is doing, but please feel free to contact me if you have any questions. I still need to sift through my code for inference and model fitting (see 2019 VSS poster) but hope to update you on that soon.

Best,
Marieke