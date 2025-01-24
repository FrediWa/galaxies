#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include<iostream>

#define GALAXIES_LENGTH 100 * 1000
#define HISTOGRAM_SIZE 90*4 // No galaxies in the data are more than 90deg apart.
#define PI 3.141592653589793238462643383

#define ARCMIN2RAD(arcmin) arcmin*PI / (60 * 180)
#define RAD2DEG(rad) rad*180 / PI

#define DD 0
#define DR 1
#define RR 2

typedef struct {
    float declination;
    float right_ascension;
} EquatorialPoint;

typedef struct
{
    int bins[HISTOGRAM_SIZE];
} Histogram;

/* ------ GPU functions ------ */
__device__ __forceinline__ float calculate_angle(EquatorialPoint* p1, EquatorialPoint* p2)
{
    float d1, d2, a1, a2;
    d1 = ARCMIN2RAD(p1->declination);
    a1 = ARCMIN2RAD(p1->right_ascension);

    d2 = ARCMIN2RAD(p2->declination);
    a2 = ARCMIN2RAD(p2->right_ascension);
    
    // Angle formula directly copied from course material.
    float b1 = sinf(d1)*sinf(d2);
    float b2 = cosf(d1)*cosf(d2)*cosf(a1 - a2);
    float b3 = b1 + b2;
    b3 = fmaxf(-1.0f, fmin(1.0f, b3)); // Clamp to -1 - 1
    float angle = acosf(b3);
    return angle;
}
/* ------ End GPU functions ------ */


__global__ void compute_angles(EquatorialPoint* galaxies1, EquatorialPoint* galaxies2, Histogram* hist, unsigned long long int* angles_calculated)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    // Get index of thread and check bound.
    if (idx >= GALAXIES_LENGTH) 
    {
        // printf("outside range\n");
        return;
    }
    // Calculate the angle between the thread's own galaxy and every other.
    for (int i = 0; i < GALAXIES_LENGTH; i++)
    {
        if (i >= GALAXIES_LENGTH) {
            printf("Thread %d: Index i out of bounds! i = %d\n", idx, i);
            return;
        }
        float angle = calculate_angle(&galaxies1[idx], &galaxies2[i]);
        angle = RAD2DEG(angle);

        // x != x if it's NaN.
        if (angle != angle) 
        {
            printf("NaN\n");
        }

        // Put into bins.
        int bin_index = (int)(angle * 4);
        
        // hist = Histogram* hist->bins
        atomicAdd(hist->bins + bin_index, 1);

        atomicAdd(angles_calculated, 1);
    }

}

/* ------ Host functions ------ */
int get_galaxy_data(EquatorialPoint* galaxy_list, int N, const char* pathname)
{
    // Open file.
    FILE* file = fopen(pathname, "r");
    
    if(file == NULL)
    {
        return 1;
    }
    // Read data from file.
    float dec = 0.0;
    float ra = 0.0;
    for (int i = 0; i < N; i++)
    {
        if(fscanf(file, "%f %f", &ra, &dec) != 2)
        {
            printf("Error reading file");
            return 1;
        }
        
        galaxy_list[i].declination = dec;
        galaxy_list[i].right_ascension = ra;
    }
    fclose(file);

    return 0;
}
/* ------ End Host functions ------ */

// Main
int main(void)
{
    printf("Entered main\n");
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    printf("Max threads per block: %d\n", prop.maxThreadsPerBlock);
    printf("Max grid size: %d\n", prop.maxGridSize[0]);

    cudaError_t      cuda_err;
    EquatorialPoint* h_galaxies_real;
    EquatorialPoint* h_galaxies_synthetic;

    EquatorialPoint* d_galaxies1;
    EquatorialPoint* d_galaxies2;

    Histogram h_hist_dd, h_hist_dr, h_hist_rr;
    Histogram* d_output;

    unsigned long long int angles_calculated = 0;
    unsigned long long int* d_angles_calculated;
    cudaMalloc(&d_angles_calculated, sizeof(unsigned long long int));
    // long int total_dd_computations, total_rr_computations, total_dr_computations;

    // Init histograms.
    memset(h_hist_dd.bins, 0, HISTOGRAM_SIZE);
    memset(h_hist_dr.bins, 0, HISTOGRAM_SIZE);
    memset(h_hist_rr.bins, 0, HISTOGRAM_SIZE);
    cudaMalloc(&d_output, sizeof(Histogram));
    cudaMemset(d_output->bins, 0, HISTOGRAM_SIZE);

    // Fetch real and synthetic data.
    h_galaxies_real = (EquatorialPoint*) malloc(GALAXIES_LENGTH * sizeof(EquatorialPoint));
    h_galaxies_synthetic = (EquatorialPoint*) malloc(GALAXIES_LENGTH * sizeof(EquatorialPoint));
    if(h_galaxies_real == NULL || h_galaxies_synthetic == NULL) return 1;

    (void) get_galaxy_data(h_galaxies_real, GALAXIES_LENGTH, "./data/data_100k_arcmin.dat");
    (void) get_galaxy_data(h_galaxies_synthetic, GALAXIES_LENGTH, "./data/rand_100k_arcmin.dat");
    printf("First galaxy of real data: %f %f\n", h_galaxies_real[0].right_ascension, h_galaxies_real[0].declination);
    // Malloc 2 galaxy lists for device.
    cuda_err = cudaMalloc(&d_galaxies1, GALAXIES_LENGTH * sizeof(EquatorialPoint));
    if(cuda_err != cudaSuccess) return 1;
    cuda_err = cudaMalloc(&d_galaxies2, GALAXIES_LENGTH * sizeof(EquatorialPoint));
    if(cuda_err != cudaSuccess) return 1;
    
    // Setup kernel.
    dim3 blockDim(1024, 1);
    dim3 gridDim((GALAXIES_LENGTH + blockDim.x - 1) / blockDim.x);

    // Setup computations and launch kernel.
    Histogram* hist;
    EquatorialPoint* h_galaxies1;
    EquatorialPoint* h_galaxies2;
    for (int hist_i = 0; hist_i < 3; hist_i++)
    {
        // Chose input galaxies and output histogram.
        switch(hist_i) {
            case DD: 
                h_galaxies1 = h_galaxies_real;
                h_galaxies2 = h_galaxies_real;
                hist = &h_hist_dd;
                break;
            case DR:
                h_galaxies1 = h_galaxies_synthetic;
                h_galaxies2 = h_galaxies_real;
                hist = &h_hist_dr;
                break;
            case RR:
                h_galaxies1 = h_galaxies_synthetic;
                h_galaxies2 = h_galaxies_synthetic;
                hist = &h_hist_rr;
                break;
        }

        printf("First h_galaxies1 declination %f\n", h_galaxies1[0].declination);
        // Copy data host -> device.
        cuda_err = cudaMemcpy(d_galaxies1, h_galaxies1, GALAXIES_LENGTH * sizeof(EquatorialPoint), cudaMemcpyHostToDevice);
        if(cuda_err != cudaSuccess)
            std::cout << "Error in memcpy" << cudaGetErrorString(cuda_err) << std::endl;
        cuda_err = cudaMemcpy(d_galaxies2, h_galaxies2, GALAXIES_LENGTH * sizeof(EquatorialPoint), cudaMemcpyHostToDevice);

        // Compute the three histograms.
        printf("Computing for %c%c\n", (hist_i < RR) ? 'D' : 'R', (hist_i > DD) ? 'R' : 'D'); // I'm both proud and disgusted by this.
        compute_angles<<<gridDim, blockDim>>>(d_galaxies1, d_galaxies2, d_output, d_angles_calculated);
        
        // Wait for a single histogram to be computed.
        cudaDeviceSynchronize();
        
        // Copy data device -> host.
        cuda_err = cudaMemcpy(hist, d_output, sizeof(Histogram), cudaMemcpyDeviceToHost);
        if(cuda_err != cudaSuccess)
            std::cout << "Error in memcpy " << cudaGetErrorString(cuda_err) << std::endl;
        cuda_err = cudaMemcpy(&angles_calculated, d_angles_calculated, sizeof(unsigned long long int), cudaMemcpyDeviceToHost);
        
        cudaMemset(d_angles_calculated, 0, sizeof(unsigned long long int));

        printf("Angles calculated %ld\n", angles_calculated);
    }

    printf("DD first value: %d\n", h_hist_dd.bins[0]);
    printf("DR first value: %d\n", h_hist_dr.bins[0]);
    printf("RR first value: %d\n", h_hist_rr.bins[0]);

    // Clean up.
    free(h_galaxies_real);
    free(h_galaxies_synthetic);
    cudaFree(d_output);
    cudaFree(d_galaxies1);
    cudaFree(d_galaxies2);
}