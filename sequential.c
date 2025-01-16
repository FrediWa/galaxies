#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEBUG 0
#define VERBOSE 1
#define VDEBUG (DEBUG && VERBOSE)

#define GALAXIES_LENGTH 1000
#define HISTOGRAM_SIZE 721
#define PI 3.1415926535

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
void save_histogram(Histogram* hist, const char* filename, const char* delimeter);
int main(void) 
{
    EquatorialPoint* real_galaxies;
    EquatorialPoint* synthethic_galaxies;
    EquatorialPoint* mock_galaxies;
    Histogram DD, DR, RR;
    init_histogram(&DD);
    init_histogram(&DR);
    init_histogram(&RR);

    // real_galaxies = malloc(GALAXIES_LENGTH * sizeof(EquatorialPoint));
    // synthethic_galaxies = malloc(GALAXIES_LENGTH * sizeof(EquatorialPoint));
    mock_galaxies = malloc(1000 * sizeof(EquatorialPoint));

    // (void) get_galaxy_data(real_galaxies, GALAXIES_LENGTH, "./data/data_100k_arcmin.dat");
    // (void) get_galaxy_data(synthethic_galaxies, GALAXIES_LENGTH, "./data/rand_100k_arcmin.dat");
    (void) get_galaxy_data(mock_galaxies, GALAXIES_LENGTH, "./data/smaller.dat");


    if(DEBUG) printf("First galaxy %f %f\n", mock_galaxies[0].declination, mock_galaxies[0].right_ascension);

    // calculate_angles(real_galaxies, synthethic_galaxies, GALAXIES_LENGTH, &DR);
    calculate_angles(mock_galaxies, mock_galaxies, GALAXIES_LENGTH, &RR);

    for (int i = 0; i < HISTOGRAM_SIZE; i++) 
    {
        if(RR.bins[i] > 0)
            printf("%d %d\n ", i, RR.bins[i]);
    }
    if(DEBUG) printf("First galaxy %f %f\n", mock_galaxies[0].declination, mock_galaxies[0].right_ascension);

    save_histogram(&RR, "output.csv", ";");

    free(mock_galaxies);
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

        if(VDEBUG) printf("Declination %f\n", dec);
    }
    fclose(file);
}

float calculate_angle(EquatorialPoint* p1, EquatorialPoint* p2, int n)
{
    float d1, d2, a1, a2;
    d1 = ARCMIN2RAD(p1->declination);
    d2 = ARCMIN2RAD(p2->declination);
    a1 = ARCMIN2RAD(p1->right_ascension);
    a2 = ARCMIN2RAD(p2->right_ascension);
    
    // Angle formula directly copied from course material.
    float b1 = sinf(d1)*sinf(d2);
    float b2 = cosf(d1)*cosf(d2)*cosf(a1 - a2);
    float b3 = b1 + b2;
    b3 = fmaxf(-1.0f, fminf(1.0f, b3));
    float angle = acosf(b3);
    // printf("Galaxy pair %f %f\n", d2, a2);
    if (VDEBUG) 
    {
        printf("%f %f %f %f\n", d1, d2, a1, a2);
        printf("b1 %f\n", b1);
        printf("b2 %f\n", b2);
        printf("b3 %f\n", b3);
        printf("acosf(1.0) = %f\n", acosf(1.0f));
        printf("acosf(b3) = %f\n", acosf(b3));
        printf("angle %f\n", angle);
    }
    return angle;
}

void calculate_angles(EquatorialPoint* galaxy_list1, EquatorialPoint* galaxy_list2, int N, Histogram* hist)
{
    for (int i = 0; i < N; i++)
    {
        printf("\rCalculating for galaxy %d", i); fflush(stdout);
        for (int j = 0; j < N; j++)
        {
            // printf("Iteration %d\n", j);
            float angle = calculate_angle(&galaxy_list1[i], &galaxy_list2[j], i);
            // printf("Galaxy pair %f %f\n", galaxy_list1[i].declination, galaxy_list1[i].right_ascension);
            // printf("Galaxy pair %f %f\n", galaxy_list2[j].declination, galaxy_list2[j].right_ascension);
            angle = RAD2DEG(angle);
            // printf("angle %f\n", angle);
            if (angle != angle) 
            {
                printf("NaN\n");
            }
            int bin_index = (int)round(angle * 4);
            // printf("What %d\n", bin_index);
            hist->bins[bin_index] += 1;
        }
    }
    printf("\nDone calculating galaxies\n");
}

void init_histogram(Histogram* hist)
{
    // Apparently gcc has undefined behaviour for initialization of struct member arrays.
    for (int i = 0; i < HISTOGRAM_SIZE; i++) 
        hist->bins[i] = 0;
}   

void save_histogram(Histogram* hist, const char* filename, const char* delimeter)
{
    FILE *file = fopen(filename, "w");

    for (int i = 0; i < HISTOGRAM_SIZE; i++) {
        printf("%d, ", hist->bins[i]);
        fprintf(file, "%d\n", hist->bins[i]); 
         
    }
}