## Understanding the problem
The first step was understanding the problem at all, what to do. 

## Data
The data is 2 lists of galaxies, one measured and one synthetic. Both input lists contain ordered pairs, right ascension and declination respectively, both in arc minutes. Arc minutes may be converted to radians with `(1/60)*pi/180`. 

It's worth noting that the data only has angles, there's no information about distance (or magnitude if we think about vectors).

# Task 
The task involves creating 3 histograms DD(real-real), DR(synthetic-real) and RR(synthetic-synthetic). Each histogram should cover 0 to 180 degrees, in bins of 0.25 degrees, for a total of 720 bins.

As per the slides, the thing that needs to be calculated is an angle $\theta_{nm} = arccos[sin(d_n)sin(d_m) + cos(d_n)cos(d_m)cos(a_n-a_m)]$. This angle will be part of the output data, visualized using a histogram. 

Optional: visualize the histogram
Mandatory: calculate $\omega_i$