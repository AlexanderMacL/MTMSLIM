# MTMSLIM

A series of MATLAB scripts to do better, more customisable, and more informed SLIM analysis.

## SLIM

The Spacer Layer Imaging Method (SLIM) is a technique to map the thickness of a thin film using white light interferometry. The last step in the process is to calculate film thicknesses based on the pixel colours in a photograph of the interference fringes. The colour matching process is simple to conceptualise but highly error-prone, and frequently produces spurious results.

This repository contains the tools to perform a more advanced colour matching, involving a context-aware colour matching which prevents aliases, a calibration curve refit step which minimises colour intensity shift, and incorporates a fast circle-finding algorithm.

The code is at present written in MATLAB and is rather dense, work to make it more user-friendly is in progress (Apr 2024) and updates will be committed to this repository as it progresses.

## How to use this code

There are four scripts of interest to the end user:
- `MTMSLIMexample.m ` runs out-the box, prompting for the user to select a `.mtmd` file, and analyses all of the SLIM images in the file, plotting the results. The script contains an example call to `MTMSLIM.m` and should be customised by the user for their purposes.
- `MTMSLIM.m` contains a function which extracts the image names from the `.mtmd` file, and is a very thin wrapper for `SLIMwrap.m`.
- `SLIMwrap.m` is a function whose arguments include the names and locations of the SLIM images for analysis, the calibration to use, and several other configuration parameters. The function returns the results of the analysis. The analysis involves several steps, including:
    - Image loading
    - A transformed-sinusoidal fit to the calibration curve (optional)
    - Circle-finding using a convolution-based algorithm to find symmetries in the image (optional)
    - An iterative refitting of this curve to the colour data from the images (optional)
    - Context-aware colour-matching using an expanding-box traverse pattern
    - Calculation and mapping of film thickness, distance and nearest calibration colours
- `SLIMwrapexample.m` which is an example usage of `SLIMwrap.m` which can be used without a `.mtmd` file, for example to analyse a specific subset of images

## Things to bear in mind

- This code is experimental
- Feedback is much appreciated (at dam216@ic.ac.uk)
- The subset of images you choose affects the results if you are using calibration curve refitting. This is because the curve is fitted to a subset of pixels within the images and so adding or removing images changes the shape of the point cloud and hence the fitted calibration curve.
- The SLIMwrap algorithm is very sensitive to the first 5 pixels in the expanding box traverse pattern, i.e. right at the centre of the image. If your image appears to have whole quadrants which are aliased, try increasing the averaging block size, or turning circle fitting off to allow manual control of the circle centre so you can avoid colour noise for these pixels.

