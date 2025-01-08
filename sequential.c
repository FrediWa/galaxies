#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEBUG 1
#define VERBOSE 0
#define VDEBUG (DEBUG && VERBOSE)

#define GALAXIES_LENGTH 100000
#define HISTOGRAM_SIZE 720
#define PI 3.1415926535

#define ARCMIN2RAD(arcmin) arcmin*PI / (60 * 180)
#define RAD2ARCMIN(rad) rad*(60*180) / PI

typedef struct {
    float declination;
    float right_ascension;
} EquatorialPoint;

typedef struct
{
    float bins[HISTOGRAM_SIZE];
} Histogram;

int get_galaxy_data(EquatorialPoint* galaxy_list, int N, const char* pathname);
float calculate_angle(EquatorialPoint* p1, EquatorialPoint* p2);
void calculate_angles(EquatorialPoint* galaxy_list1, EquatorialPoint* galaxy_list2, int N, Histogram* hist);
void init_histogram(Histogram* hist);

int main(void) 
{
    EquatorialPoint* real_galaxies;
    EquatorialPoint* synthethic_galaxies;
    Histogram DD, DR, RR;
    init_histogram(&DD);
    init_histogram(&DR);
    init_histogram(&RR);

    real_galaxies = malloc(GALAXIES_LENGTH * sizeof(EquatorialPoint));
    synthethic_galaxies = malloc(GALAXIES_LENGTH * sizeof(EquatorialPoint));

    (void) get_galaxy_data(real_galaxies, GALAXIES_LENGTH, "./data/data_100k_arcmin.dat");
    (void) get_galaxy_data(synthethic_galaxies, GALAXIES_LENGTH, "./data/rand_100k_arcmin.dat");
    
    if(DEBUG) printf("First galaxy %f %f\n", real_galaxies[0].declination, real_galaxies[0].right_ascension);

    free(real_galaxies);
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
        if(fscanf(file, "%f %f", &dec, &ra) != 2)
        {
            if(DEBUG) printf("Error reading file");
            return 1;
        }
        
        galaxy_list[i].declination = dec;
        galaxy_list[i].right_ascension = ra;

        if(VDEBUG) printf("%f\n", galaxy_list[i].declination);
    }
    fclose(file);
}

float calculate_angle(EquatorialPoint* p1, EquatorialPoint* p2)
{
    float d1, d2, a1, a2;
    d1 = p1->declination;
    d2 = p2->declination;
    a1 = p1->right_ascension;
    a1 = p2->right_ascension;

    // Angle formula directly copied from course material.
    float angle = acosf(sinf(d1)*sinf(d2) + cosf(d1)*cosf(d2)*cosf(a1 - a2));
}

void calculate_angles(EquatorialPoint* galaxy_list1, EquatorialPoint* galaxy_list2, int N, Histogram* hist)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            float angle = calculate_angle(&galaxy_list1[i], &galaxy_list2[j]);
        }
    }
}

void init_histogram(Histogram* hist)
{
    // Apparently gcc has undefined behaviour for initialization of struct member arrays.
    for (int i = 0; i < HISTOGRAM_SIZE; i++) 
        hist->bins[i] = 0;
}   


/*
0 | 0.25 | 0.5 | 0.75 | 1 | 1.25
0.37 -> 1
1.1-> 5


*/