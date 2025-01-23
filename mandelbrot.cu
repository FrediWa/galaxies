#include <iostream>
#include <cuda_runtime.h>
#include <cuComplex.h>

#define ITERATIONS 50
#define RESOLUTION 10000
#define TOTALSIZE RESOLUTION*RESOLUTION 

__device__ __forceinline__ cuComplex csquare (cuComplex z)
{
    cuComplex res;
    res.x = (z.x*z.x - z.y*z.y);
    res.y = (2*z.x*z.y);
    return res;
}

__global__ void mandelbrot(char* output) {
        // Perform computation here
    int xidx = threadIdx.x + blockIdx.x * blockDim.x;
    int yidx = threadIdx.y + blockIdx.y * blockDim.y;
    
    if (xidx >= RESOLUTION || yidx >= RESOLUTION) {
        return;
    }
    // Create the starting complex number.
    cuComplex c = make_cuComplex(xidx, yidx);

    // Perform Mandelbrot computation.
    cuComplex z;
    for (int i = 0; i < ITERATIONS; i++)
    {
        z = csquare(z);
        z.x += c.x;
        z.y += c.y;
    }

    // Assign value based on magnitude
    float magnitude = cuCabsf(z);
    char single_result;
    if (magnitude < 2) 
        single_result = 'a';
    else if (magnitude < 5) 
        single_result = 'e';
    else if (magnitude < 20) 
        single_result = 'c';
    else if (magnitude < 100) 
        single_result = 'd';
    else 
        single_result = 'e';

    int idx = yidx * RESOLUTION + xidx;
    
    output[idx] = single_result;
}

int main(void)
{
    std::cout << "Entered main" << std::endl;
    dim3 blockDim(32, 32); 
    dim3 gridDim((RESOLUTION + blockDim.x - 1) / blockDim.x, 
                (RESOLUTION + blockDim.y - 1) / blockDim.y); 

    // Allocate temp memory on device.
    char* d_output;
    cudaMalloc(&d_output, TOTALSIZE * sizeof(char));
    // Allocate output memory on host.
    char* output = (char*) malloc(TOTALSIZE * sizeof(char));
    if(output == NULL)
    {
        std::cout << "error mallocing" << std::endl;
        return -1;
    }

    // Compute the mandelbrot set points.
    mandelbrot<<<gridDim, blockDim>>>(output);
    cudaError_t err = cudaMemcpy(output, d_output, TOTALSIZE * sizeof(char), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        std::cout << "cudaMemcpy failed: " << cudaGetErrorString(err) << std::endl;
        return -1;
    }

    std::cout << "[" << output[0] << "]" << std::endl;
    // Write output to file.
    // FILE* f = fopen("output.txt", "w");
    // for(int i = 0; i < TOTALSIZE; i++)
    // {
    //     std::cout << 'a';
    // }
    // fclose(f);
    std::cout << "came here" << std::endl;

    // Cleanup.
    free(output);
    cudaFree(d_output);
}