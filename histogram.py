import csv
import matplotlib.pyplot as plt
import np

bins = np.arange(0, 180.25, 0.25)

with open('output.csv', newline='') as f:
    reader = csv.reader(f, delimiter=';')
    data = [int(row[0]) for row in reader] 

plt.bar(bins, data, width=0.25, edgecolor='black', color='skyblue')
plt.show()