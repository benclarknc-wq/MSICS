# MSICS
Multi Speed Image Correlation Spectroscopy (MSICS)

Problem: Deadtime in RICS analysis lead to measurement gap in image analysis.

Solution: Use Bidirectional Scanning to have continuous collection

Result: New image analysis tool

MSICS is a image correlation tool aims to accurately measure the diffusion dynamics of molecules from bidirectional mircoscope images of flourescent molecules. This is accomplished by conducting Spatio-Temporal Image Correlation Spectroscopy (STICS), which is a class of image analysis that compares photon counts throughout the pixels of an image based on the time lag between collected pixels and the spatial displacement between pixels. Other STICS methods exist looking between individual image frames (ICS probes on the order of ms), raster images (RICS probes on the order of 10 us), or immobile time traces (FCS probes on the order of 10 ns). Each method can probe different time scales dependent upon the geometery of the image collected. MSICS aims to probe time scales on the us scale a happy medium between RICS and FCS. By conducting this analysis one can gain diffusion information, the diffusion coefficeint of molecules in a solution, as well as the molecule concentration in the solution. This analysis was conducted on a prepared solution of 10uM streptavidin-quantum dot complex in an oxygen scavenging imaging buffer that emit light at the wavelength of 580nm (10uM QD 580). A 100x100 pixel image of 600 frames with a pixel size of 30 nm were collected in this example, the solution was excited with a 485 nm laser with a confocal radius of 353 um and aspect ration of ~3.333. 30 nm was chosen as the pixel size to over sample the excitation confocal volume.

<img width="600" height="300" alt="image" src="https://github.com/user-attachments/assets/bafbe56c-6d61-4a72-ada5-cd214afceb63" />

The figure on the left depicts a raster scan the black arrows indicate times when photon collection is off, while the figure on the right depicts a bidirectional scan. The bidirectional scan allows for continuous collection of photons with out collection gaps. Once an image stack is collected an auto correlation is conducted on each image frame with respect to $\xi$ pixel collumn shifts and $\psi$ pixel row shifts and is then averaged across the image stack as shown below.

<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/353e0d7d-33a8-44b9-a7c4-026216589523" />

A unique trait of bidirectional images is distortion due to a phenomonon known as hysteresis, which in this case refers to lag in motion of the scanning laser due to charge build up on the piezo mirror which drives the laser. The lag causes noticeable streaking in output images which we will correct in this analysis as well using a sinesoidal model to capture the lagging motion. An example of both the hysteresis behavior and correction are as follows:

<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/3f2bc6fd-6ae4-433f-93ae-1288be7523d0" />
<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/095535ce-12de-4b3e-8971-1e169e953fcc" />

This image was taken from fixed flourescent beads, the beads should form a single circular reponse in the resuling image however there are noticable streaks which form between even and odd rows. Multiple images where taken of these immobile beads through the imaging field of a bidirectional scan and transposed on top of one another to find correction parameters.

Once image frames are correlated and spatially corrected for hysteresis we produce an ROC curve which depicts correlation strength relative to timelag over many different $\xi$ and $\psi$ shifts. The following is an example ROC curve which considers purely $\xi$ shifts and $\psi$ shifts: 

<img width="540" height="380" alt="image" src="https://github.com/user-attachments/assets/7d0d64f8-01de-4c84-bb6c-3bcb1cebc4aa" />

This ROC can be fit by the following equations which desribe the natural diffusion of the particles and scanning behavior of the laser:

$G(\tau(\xi,\psi),\Delta r) = G_D(\tau(\xi,\psi))S(\Delta r)$

Broken down further into their more elementary formulas we find:

$G_D(\tau) = \frac{\gamma}{\bar{N}}\left(1+\frac{4D|\tau(\xi,\psi)|}{\omega_0^2}\right)^{-1}\left(1+\frac{4D|\tau(\xi,\psi)|}{\kappa^2\omega_0^2}\right)^{-1/2}$

$S(\Delta r) = \exp{\left[-\frac{\left(\frac{|\xi|\delta x}{\omega_0}\right)^2+\left(\frac{|\psi|\delta y}{\omega_0}\right)^2}{1+\frac{4D|\tau(\xi,\psi)|}{\omega_0^2}}\right]} $

Here, $\omega_0$ is the lateral radius of the scanning laser, $\kappa$ is the aspcet ratio of the confocal volume. The variable $\gamma \approx 0.3535$ is a correction factor for uneven illunimation throughout a 3D gaussian confocal volume. The pixel size is described by $\delta x,\delta y$. As timelag (and spatial displacement) increase the correlation of photon signal should decrease, that is to say a diffusing particle will move and the light signal we detect from that particle will fall over time. The two fit parameters are the target of the analysis: D the diffusion coefficient and N which is related to the concentration of the solution. An example of this fitting process is shown below:

<img width="900" height="600" alt="image" src="https://github.com/user-attachments/assets/c1580d4f-e3cd-4427-bc2f-50eabc176b63" />

This ROC curve has both the fit curves with and without hysteresis correction. The hysteresis correction helps explain some of the characteristic dips which occur thoughout the ROC curve. In this case the fit D and N values for the streptavidin-QD complexes where 21.4 $\mu m^2/s$ and 0.7. Which matches the measured diffusion coefficient via FCS (22.0 $\mu m^2/s$) and from theoretical diffusion coefficients of QD580 in ddH20 (24 $\mu m^2/s$).

This new correlation method expands the probable rages for ICS methods being able to reach lag times which RICS cannot access as well as the spatial information which FCS cannot access. There is still work to be done on some key behaviors of this method in particular when scanning speeds begin to reach fast speeds (<1.7 us) and slow speeds (>20 us). In these cases other methods may be better served to test the diffusion coeffiecent of these diffusing particles. But over all this method is robust and matches theory as a new technique for dynamic image analysis.



