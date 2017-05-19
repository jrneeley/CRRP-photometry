# CRRP-photometry

Under construction...


Pipeline to perform photometry on IRAC images of crowded fields, using an input catalog derived from
higher resolution data. 


Pipeline.py

1. Converts images from flux units (MJy/sr/sec) to DN
2. Performs initial DAOPHOT photometry (find and phot) on all images
3. Copies master PSF to all images - you must create the PSF ahead of time
4. Runs ALLSTAR
5. Runs DAOMATCH between all of your individual frames
6. Refines the transformations using DAOMASTER
7. Identifies the appropriate x and y values in the optical catalog for each of your fields.

