# AMBER V1.0
**Artefactual Multiphoton Bundle Effect Removal**
## Description
This release is a collection of functions and an example of a main program that is intended for the use of removing an artifact caused by performing multiphoton imaging through a coherent imaging fiber bundle. Some of the code may need to be modified to be compatible with specific file types for the image stacks. Also required is a "Flat" image, which is a multiphoton image of a flat fluorescent sample through the fiber-bundle acquired under the same conditions as the actual image data stacks to be fixed. The individual image stacks may be time series or Z-stacks and may be a field-of-view smaller than the full surface of the fiber-bundle. A simple rigid registration code allows the centroid locations identified in the flat image to be registered onto the image data.
## List of functions:
- centerPadCrop: Performs a simple central crop
- filterValSeries: Temporal/axial filter on each core independently
- getCentroidValues: Extracts the values of the cores at the identified centroids
- getCorrectorValues: Acquires the correction factors for each core from the flat image
- getFiberCentroids: Finds the centroid pixel coordinates for each fiber core
- getParameter: Extract specific data from a *.txt file containing imaging metadata
- gridFiberCores: Interpolates the core centroids into a uniform pixel grid
- makeFiberImage: Uses the centroid locations and values to recreate the fiber images
- rigidAlignFiber: Performs a rigid registration of the flat fiber centroids to the actual image field
- writeFiberImages: Writes the processed image files to *.tif with metadata (using Bioformats toolbox)

Author: Baris N. Ozbay, University of Colorado Denver, Department of Bioengineering

https://zenodo.org/badge/DOI/10.5281/zenodo.1247507.svg
