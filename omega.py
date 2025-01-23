
histograms = {"dd": list(), "dr": list(), "rr": list(), }
DD = list()
DR = list()
RR = list()

for hist in ["dd", "dr", "rr"]:
    f = open(f"{hist}.csv", "r")
    for line in f:
        match hist:
            case "dd": DD.append(int(line))
            case "dr": DR.append(int(line))
            case "rr": RR.append(int(line))
    f.close()

omega = list()
totalDD = 0
totalDR = 0
totalRR = 0

for i in range(0, 360):
    omega.append((DD[i] - 2 * DR[i] + RR[i]) / RR[i])
    totalDD += DD[i]
    totalDR += DR[i]
    totalRR += RR[i]

print(omega)
print("totalDR", totalDR)
print("totalDD", totalDD)
print("totalRR", totalRR)
