# MSICS
Multi Speed Image Correlation Spectroscopy

This image correlation tool aims to accurately measure the diffusion dynamics of molecules from bidirectional mircoscope images of flourescent molecules. This is accomplished by conducting Spatio-Temporal Image Correlation Spectroscopy (STICS), which is a class of image analysis that compares photon counts throughout the pixels of an image based on the time lag between collected pixels and the spatial displacement between pixels. Other STICS methods exist looking between individual image frames (ICS probes one the order of ms), raster images (RICS probes on the order of 10 us), or immobile time traces (FCS probes on the order of 10 ns). Each method can probe different time scales dependent upon the geometery of the image collected. MSICS aims to probe time scales on the us scale a happy medium between RICS and FCS. By conducting this analysis one can gain diffusion information, the diffusion coefficeint of molecules in a solution, as well as the molecule concentration in the solution.



This analysis was conducted on a prepared solution of 10uM quantum dots in an oxygen scavenging imaging buffer that emit light at the wavelength of 580nm (10uM QD 580). A 100x100 pixel image of 600 frames with a pixel size of 30 nm were collected in this example, the solution was excited with a 485 nm laser with a confocal radius of 353 um and aspect ration of ~3.333. 30 nm was chosen as the pixel size to over sample the excitation confocal volume.

A unique trait of bidirectional images is distortion due to a phenomonon known as hysteresis, which in this case refers to lag in motion of the scanning laser due to charge build up on the piezo mirror which drives the laser. The lag causes noticeable streaking in output images which we will correct in this analysis as well using a sinesoidal model to capture the lagging motion.

The output of the analysis is an ROC curve which depicts correlation strength relative to timelag. This correlation follows from theory and can be described by the following equation:

$G(\tau,\Delta r) = G_D(\tau)S(\Delta r)$

Broken down further into their more elementary formulas we find:

$G_D(\tau) = \frac{\gamma}{\bar{N}}\left(1+\frac{4D|\tau(\xi,\psi)|}{\omega_0^2}\right)^{-1}\left(1+\frac{4D|\tau(\xi,\psi)|}{\kappa^2\omega_0^2}\right)^{-1/2}$

$S(\Delta r) = \exp{\left[-\frac{\left(\frac{|\xi|\delta x}{\omega_0}\right)^2+\left(\frac{|\psi|\delta y}{\omega_0}\right)^2}{1+\frac{4D|\tau(\xi,\psi)|}{\omega_0^2}}\right]} $

As timelag (and spatial displacement) increase the correlation of photon signal should decrease, that is to say a diffusing particle will move and the light signal we detect from that particle will fall over time. We can fit this ROC curve with respect to the natural diffusion of particles G(tau) and relative to the scanning displacement of the laser S(r). The two fit parameters are the target of the analysis: D the diffusion coefficient and C the concentration of the solution.

