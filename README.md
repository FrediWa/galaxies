# Galaxy calculator

## Running
For sequntial run `make` followed by `a.out`.

Parallel is compiled with `nvcc -arch=sm_70 galaxy.cu -o galaxy -lm` and run with `srun -p gpu -n 1 -t 10:00 --mem=1G -e err.txt galaxy`. 

## Introduction
The problem is to design and implement a program that calculates omega values for a distribution of galaxies. With 100000 random galaxies and 100000 real galaxies, times 3 histograms, the number of computations to be run are 30 billion. 

## Computations for each galaxy
For each galaxy, an angle between it and every other galaxy is to be computed. For each galaxy, read every other galaxy's coordinates and compute the angle. The angle formula is derived in the material and fairly easy to implement. This angle, in radians should be turned into degrees and then put into quarter-degree histogram bins.

The data is such that no galaxy has an angle of 90 degrees of greater. This means that a total of 360 bins is needed for the histogram.

## Main
There seems to be a common setup for the main and consists of the following steps:
- Use cudaMemcpy to copy data from host to device
- Perform computations
- Use cudaMemcpy to copy data from device to host

Other things that are common to cuda programs: allocate the data on both the host and device with malloc and cudaMalloc respectively (and freeing them), setting up the grid and block dimensions. When I computed the Mandelbrot set as a warmup exercise, I used a 2D block and grid which made indexing in x and y dimensions possible. For the galaxies, it makes probably more sense to only have one dimension as the data is also a single dimensions (even though every data point has 2 components). Int the case of 2D, I used 32x32, but for the galaxies I instead use 1024x1. 

## Kernel function
The kernel functions also seems to have a fairly standard template: calculate an id based on the block and thread index then do something and perform I/O based on the index. Then since the range of the computations are seldom a multiple of the block size (in any dimension), a bounds check is added to make sure that only threads with indexes within the range are actually doing I/O. If they did computations, that's at worse reduntant, but I/O outside bounds is worse. A Cuda refresher article on developer.nvidia.com has an example of matrix multiplication. The threads are assigned IDs in two dimensions, then every thread is responsible for calculating a single entry in the output matrix. The bound is checked by making sure that a threads x-dimension ID is less than the width of the matrix and same is done for the y-dimension ID.

The neat thing with the computations is that almost everything happens within a thread. The computations involves reading 2 equatorial point, computing the angle between the them and writing to a histogram. With the Mandelbrot set, every thread was writing to essentially it's own region of memory (really just a region based on it's own position within the block and the grid), but the galaxy calculations need to write to a histogram, which is shared among the threads. The histogram bins thus need to be incremented atomically.

## Sequential implemention
The sequential implementation was fairly straightforward, but I had a lot of issues in getting it to work. I got together something that should have worked in an hour, but then it didn't work and after a long time, it was pointed out that my histogram index calculations are wrong. The omega calculations were also wrong, but I'm guessing I would've found them had I atleast gotten the right histograms. Oh well.

## Mandelbrot set
In the interim of waiting for help, I wanted to try to write another Cuda program. I asked ChatGPT for help and it gave me a bunch of stuff like matrix multiplication (which is like a GPU hello world). It then suggested the mandelbrot set which I thought would be cool. It only took me a couple of hours to get a working mandelbrot set program. As Dione is a server, my solution was to print characters instead of value. A python program would then visualize the set. 

It was a good exercise because I fell for some common-ish? pitfalls. It took me at least an hour to find out why my cudaMemcpy was failing. The reason was that I gave the kernel the host's allocated memory and not the device's allocated memory.

## Galaxy computations
After struggling with stupid mistakes for with the sequential implementation, I migrated all code to a cuda program. The warmup exercise with the mandelbrot really helped me in creating the Cuda galaxy implementation and I did have some pitfalls (probably due to trying to be smart with how I set up the computations). Like with the sequential program, I first wanted to see that the histograms contained the right data. Weirdly, they did not, even though the program was almost a direct copy from the working sequential program. The problem was that I used the same Histogram* d_output on the device, and I forgot to reset it between runs. Adding a cudaMemset to reset d_output fixed the problem and now gives more or less the right histograms. 

## Computing omega
When it came to computing omega, I had to do some restructuring of the program. Even though the omega computations (the histograms being small in size and not dependant on input size) wouldn't benefit from parallelization, they were still very easy to parallelize. First thing I realized that when I previously had 3 histograms on the host and 1 on the device, I now wanted the opposite because the 3 histograms were needed simultaneously on the device. I now allocated 3 device histograms and only one host histogram (to print that atleast the first bin has correct values). Once the angle calculations were done, I launched another kernel with the typical setup (except no host to device memcpy because the input is already on the device). The same grid/block dimensions were used and the output was memcpy'ed deviceToHost. The final output, first 5 of which are printed, match the values on moodle. 

## Conclusion
Initially, the problem was hard to understand. After a while I got confident that I understood the problem so I started writing a sequential program. I had lots of troubles with this but it really just boiled down to incorrectly calculating the bin_index for a certain galaxy pair. While thinking of this I warmed myself up to Cuda by computing the mandelbrot set. After getting that to work, I started parallelizing the galaxies. I struggled a bit with it, but the most common pitfalls were already behind me from the mandelbrot program. Following, what seems like a fairly standard template for Cuda programs, I got the computations fairly easy migrated to Cuda. 

