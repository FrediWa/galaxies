
histograms = {"dd": list(), "dr": list(), "rr": list(), }
DD = list()
DR = list()
RR = list()

for hist in ["dd", "dr", "rr"]:
    f = open(f"{hist}.csv", "r")
    for line in f:
        # print(line.rstrip())
        
        match hist:
            case "dd": DD.append(int(line))
            case "dr": DR.append(int(line))
            case "rr": RR.append(int(line))
    f.close()

print(len(DD))
# for line in histograms["dd"]:
#     print(line)

omega = list()

for i in range(0, 361):
    omega.append((DD[i] - 2 * DR[i] + RR[i]) / RR[i])

print(omega)
