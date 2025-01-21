#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEBUG 0
#define VERBOSE 1
#define VDEBUG (DEBUG && VERBOSE)

#define DD 1
#define DR 2
#define RR 3
#define TARGET RR

#define GALAXIES_LENGTH 100 * 1000
// No galaxies in the data are more than 90deg apart.
#define HISTOGRAM_SIZE 361 
#define PI 3.141592653589793238462643383

#define ARCMIN2RAD(arcmin) arcmin*PI / (60 * 180)
#define RAD2DEG(rad) rad*180 / PI

typedef struct {
    double declination;
    double right_ascension;
} EquatorialPoint;

typedef struct
{
    int bins[HISTOGRAM_SIZE];
} Histogram;

int get_galaxy_data(EquatorialPoint* galaxy_list, int N, const char* pathname);
double calculate_angle(EquatorialPoint* p1, EquatorialPoint* p2, int n);
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
            calculate_angles(synthethic_galaxies, synthethic_galaxies, GALAXIES_LENGTH, &hist);
            save_histogram(&hist, "dd.csv"); break;
        case DR:
            calculate_angles(real_galaxies, synthethic_galaxies, GALAXIES_LENGTH, &hist);
            save_histogram(&hist, "dr.csv"); break;
        case RR:
            calculate_angles(real_galaxies, real_galaxies, GALAXIES_LENGTH, &hist);
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
    double dec = 0.0;
    double ra = 0.0;
    for (int i = 0; i < N; i++)
    {
        if(fscanf(file, "%lf %lf", &ra, &dec) != 2)
        {
            if(DEBUG) printf("Error reading file");
            return 1;
        }
        
        galaxy_list[i].declination = dec;
        galaxy_list[i].right_ascension = ra;
    }
    fclose(file);
}

double calculate_angle(EquatorialPoint* p1, EquatorialPoint* p2, int n)
{
    double d1, d2, a1, a2;
    d1 = ARCMIN2RAD(p1->declination);
    a1 = ARCMIN2RAD(p1->right_ascension);

    d2 = ARCMIN2RAD(p2->declination);
    a2 = ARCMIN2RAD(p2->right_ascension);
    
    // Angle formula directly copied from course material.
    double b1 = sin(d1)*sin(d2);
    double b2 = cos(d1)*cos(d2)*cos(a1 - a2);
    double b3 = b1 + b2;
    b3 = fmax(-1.0f, fmin(1.0f, b3)); // Clamp to -1 - 1
    double angle = acos(b3);
    return angle;
}

void calculate_angles(EquatorialPoint* galaxy_list1, EquatorialPoint* galaxy_list2, int N, Histogram* hist)
{
    long int angles_calculated = 0;
    for (int i = 0; i < N; i++)
    {
        printf("\rCalculating for galaxy %d", i); fflush(stdout);
        for (int j = 0; j < N; j++)
        {
            double angle = calculate_angle(&galaxy_list1[i], &galaxy_list2[j], i);
            angle = RAD2DEG(angle);

            // x != x if it's NaN.
            if (angle != angle) 
            {
                printf("NaN\n");
            }
            int bin_index = (int)round(angle * 4);
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
    printf("Histogram total: %ld\n", acc);
    printf("Total computations should be %ld\n", ((long)GALAXIES_LENGTH*(long)GALAXIES_LENGTH));
}