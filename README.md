# Split-Bregman-ST-Total-Variation-MRI
Split Bregman spatiotemporal total variation for MRI

This repository contains a demo to use Spatiotemporal Total Variaton (ST-TV) solved using the Split Bregman formulation for cardiac cine MRI
used in the paper: 
**P Montesinos, J F P J Abascal, L Cussó, J J Vaquero, M Desco. Application of the compressed sensing technique to self-gated cardiac cine sequences in small animals. Magn Reson Med., 72(2): 369–380, 2013.** 
DOI: http://dx.doi.org/10.1002/mrm.24936

The demo uses prospective cardiac cine small-animal data to simulate an undersampling pattern based on a variable density pdf and compare Spatial TV with Spatiotemporal TV. Both methods are efficiently solved in the Fourier domain with a computational cost of three FFT at each iteration. 

![Image of Yaktocat](https://github.com/HGGM-LIM/Split-Bregman-ST-Total-Variation-MRI/blob/master/image_cardiac_cine. gif)

The repository contains the following files:

- **ProspectiveCine16fr.mat:** Prospective cardiac cine small-animal data (16 frames, healthy rat)

- **image_HighDose.gif:** Video of images for high dose that shows respiratory motion

- **Demo_SpatioTemporalTV_SplitBregman.m:** This demo loads data, apply an undersampling and compares S-TV and ST-TV

- **SpatialTV_SplitBregman.m:** Split Bregman S-TV method
- 
- **SpatioTemporalTVSBIso.m:** Split Bregman ST-TV 

Data has been saved using MATLAB R2012b.

If you need to contact the author, please do so at juanabascal78 (AT) gmail (DOT) com.
