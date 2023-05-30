# Active Surfaces: Chan-Vese Segmentation

![Alt text](/imgs/level_set.png "")

![Alt text](/imgs/lambdas.png "")

The Chan-Vese segmentation algorithm is a classic technique that allows segmenting objects without clearly defined boundaries
and made up of potentially several disconnected regions.

This is achieved by formulating the solution as the intersection of the image with a level set function.
It is this level set function that is evolved over time to minimize an energy, which involves a penalty on contour length (to make boundaries smooth)
and two data fidelity terms (that enforce homogeneity) penalizing deviation from the mean values inside and outside the segmented regions.

The algorithm is described in detail in the original paper 'Chan-Vese Segmentation' by Pascal Getreuer, [2012](https://www.ipol.im/pub/art/2012/g-cv).


## Running the code

Open a matlab console, and run `start.m`.

You will need to specify the path to the desired image, as well as the initial contour and hyperparameters.
The code will then run the segmentation algorithm and display the results.
