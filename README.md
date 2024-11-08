Tenue is a package for reducing and stacking astronomical images.

Tenue does lots of things wrong, including:

- You have to manually select a star for fine alignment.
- It does not do optimal weighting.
- You have to manually select the fraction of images you wish to reject.
- It makes one sky image for all of the images in the stack.
- It doesn't even make that sky image correctly, but just takes the median of the stack after subtracing a constant sky from each image.
- It doesn't iteratively mask sources in the sky.
