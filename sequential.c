#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEBUG 0
#define VERBOSE 1
#define VDEBUG (DEBUG && VERBOSE)

#define DD 1
#define DR 2
#define RR 3
#define TARGET DR

#define GALAXIES_LENGTH 100000
// No galaxies in the data are more than 90deg apart.
#define HISTOGRAM_SIZE 360
#define PI 3.141592653589793238462643383

#define ARCMIN2RAD(arcmin) arcmin*PI / (60 * 180)
#define RAD2DEG(rad) rad*180 / PI

typedef struct {
    float declination;
    float right_ascension;
} EquatorialPoint;

typedef struct
{
    int bins[HISTOGRAM_SIZE];
} Histogram;

int get_galaxy_data(EquatorialPoint* galaxy_list, int N, const char* pathname);
float calculate_angle(EquatorialPoint* p1, EquatorialPoint* p2, int n);
void calculate_angles(EquatorialPoint* galaxy_list1, EquatorialPoint* galaxy_list2, int N, Histogram* hist);
void init_histogram(Histogram* hist);
void save_histogram(Histogram* hist, const char* filename);
int main(void) 
{
    EquatorialPoint* real_galaxies;
    EquatorialPoint* synthethic_galaxies;

    Histogram hist;
    init_histogram(&hist);

    real_galaxies = malloc(GALAXIES_LENGTH * sizeof(EquatorialPoint));
    synthethic_galaxies = malloc(GALAXIES_LENGTH * sizeof(EquatorialPoint));

    (void) get_galaxy_data(real_galaxies, GALAXIES_LENGTH, "./data/data_100k_arcmin.dat");
    (void) get_galaxy_data(synthethic_galaxies, GALAXIES_LENGTH, "./data/rand_100k_arcmin.dat");

    switch(TARGET) 
    {
        case DD:
            calculate_angles(real_galaxies, real_galaxies, GALAXIES_LENGTH, &hist);
            save_histogram(&hist, "dd.csv"); break;
        case DR:
            calculate_angles(real_galaxies, synthethic_galaxies, GALAXIES_LENGTH, &hist);
            save_histogram(&hist, "dr.csv"); break;
        case RR:
            calculate_angles(synthethic_galaxies, synthethic_galaxies, GALAXIES_LENGTH, &hist);
            save_histogram(&hist, "rr.csv"); break;
            

    }

    printf("First galaxy %f %f\n", synthethic_galaxies[0].declination, synthethic_galaxies[0].right_ascension);

    free(real_galaxies);
    free(synthethic_galaxies);

    return(0);
}

int get_galaxy_data(EquatorialPoint* galaxy_list, int N, const char* pathname)
{
    // Open file.
    FILE* file = fopen(pathname, "r");
    if(VDEBUG) ("%s\n", pathname);
    

    if(file == NULL)
    {
        if(DEBUG) ("Error opening file");
        return 1;
    }
    // Read data from file.
    float dec = 0.0;
    float ra = 0.0;
    for (int i = 0; i < N; i++)
    {
        if(fscanf(file, "%f %f", &ra, &dec) != 2)
        {
            if(DEBUG) printf("Error reading file");
            return 1;
        }
        
        galaxy_list[i].declination = dec;
        galaxy_list[i].right_ascension = ra;
    }
    fclose(file);
}

float calculate_angle(EquatorialPoint* p1, EquatorialPoint* p2, int n)
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

void calculate_angles(EquatorialPoint* galaxy_list1, EquatorialPoint* galaxy_list2, int N, Histogram* hist)
{
    long int angles_calculated = 0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            float angle = calculate_angle(&galaxy_list1[i], &galaxy_list2[j], i);
            angle = RAD2DEG(angle);

            // x != x if it's NaN.
            if (angle != angle) 
            {
                printf("NaN\n");
            }
            int bin_index = (int)(angle * 4);
            hist->bins[bin_index]++;
            angles_calculated++;
        }
    }
    printf("\nDone calculating %ld galaxies\n", angles_calculated);
}

void init_histogram(Histogram* hist)
{
    // Apparently gcc has undefined behaviour for initialization of struct member arrays.
    for (int i = 0; i < HISTOGRAM_SIZE; i++) 
        hist->bins[i] = 0;
}   

void save_histogram(Histogram* hist, const char* filename)
{
    FILE *file = fopen(filename, "w");
    long int acc = 0;
    for (int i = 0; i < HISTOGRAM_SIZE; i++) 
    {
        fprintf(file, "%d\n", hist->bins[i]);
        acc += hist->bins[i];
    }
    fclose(file);
    
    printf("Histogram total: %ld\n", acc);
    printf("Total computations should be %ld\n", ((long)GALAXIES_LENGTH*(long)GALAXIES_LENGTH));
}