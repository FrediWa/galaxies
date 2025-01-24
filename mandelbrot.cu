#include <iostream>
#include <cuda_runtime.h>
#include <cuComplex.h>

#define ITERATIONS 30
#define RESOLUTION 2000
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
    
    int idx = yidx * RESOLUTION + xidx;
    if (idx >= TOTALSIZE) {
        return;
    }
    
    float x_val = (2.0f * xidx) / RESOLUTION - 1.0f;
    x_val -= 0.5f;
    float y_val = (2.0f * yidx) / RESOLUTION - 1.0f; 

    // printf("Thread (%d, %d): idx=%d\n", xidx, yidx, idx);

    // Create the starting complex number.
    cuComplex c = make_cuComplex(x_val, y_val);
    
    // Perform Mandelbrot computation.
    cuComplex z = make_cuComplex(0.0f, 0.0f);
    for (int i = 0; i < ITERATIONS; i++)
    {
        if (z.x > 10000 || z.y > 10000)
        {
            continue;
        }
        
        z = csquare(z); 
        z.x += c.x;
        z.y += c.y;

    }
    // printf("z = (%f + %fi)\n", z.x, z.y);

    // Assign value based on magnitude
    float magnitude = cuCabsf(z);
    // printf("mag = %f\n", magnitude);
    char single_result;
    if (magnitude < 2) 
        single_result = 'a';
    else if (magnitude < 500) 
        single_result = 'e';
    else if (magnitude < 2000) 
        single_result = 'c';
    else if (magnitude < 10000) 
        single_result = 'd';
    else 
        single_result = 'e';
    
    output[idx] = single_result;
    __syncthreads();
}

int main(void)
{

    cudaError_t err;

    dim3 blockDim(32, 32); 
    dim3 gridDim((RESOLUTION + blockDim.x - 1) / blockDim.x, 
                (RESOLUTION + blockDim.y - 1) / blockDim.y); 

    // Allocate temp memory on device.
    char* d_output;
    err = cudaMalloc(&d_output, TOTALSIZE * sizeof(char));
    if (err != cudaSuccess) {
        std::cout << "Cuda malloc failed: " << cudaGetErrorString(err) << std::endl;
        return -1;
    }
    cudaMemset(d_output, '0', TOTALSIZE * sizeof(char));

    // Allocate output memory on host.
    char* output = (char*) malloc(TOTALSIZE * sizeof(char));
    if(output == NULL)
    {
        std::cout << "Error mallocing" << std::endl;
        return -1;
    }

    // Compute the mandelbrot set points.
    mandelbrot<<<gridDim, blockDim>>>(d_output);
    cudaError_t kernelErr = cudaGetLastError();
    if (kernelErr != cudaSuccess) {
        std::cout << "Kernel launch error: " << cudaGetErrorString(kernelErr) << std::endl;
        return -1;
    }
    cudaDeviceSynchronize();

    err = cudaMemcpy(output, d_output, TOTALSIZE * sizeof(char), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        std::cout << "Cuda memcpy failed: " << cudaGetErrorString(err) << std::endl;
        return -1;
    }
    for (int i = 0; i < TOTALSIZE; i++)
    {
        printf("%c", output[i]);
    }

    // Cleanup.
    free(output);
    cudaFree(d_output);
}