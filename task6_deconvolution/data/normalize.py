# this script normalizes data

import sys

filename = sys.argv[1]

lines = []
with open(filename, "r") as fin:
    for line in fin:
        lines.append(line)

points = []
for line in lines:
    line = line.strip()
    # print(line)
    parts = line.split()
    # print(parts)
    x = parts[0]
    y = parts[1]
    x = float(x)
    if "^" not in y:
        y = float(y)
    else:
        base_exp = y.split('*')
        assert len(base_exp) == 2
        base, exp = base_exp
        assert exp.startswith("10^")
        exp = exp[3:]
        # print(base, exp)
        base = float(base)
        exp = float(exp)
        y = base * 10**exp
    # print(x, y)
    points.append((x, y))

max_y = max(map(lambda p: p[1], points))

points_new = []
for point in points:
    x, y = point
    y_new = y / max_y
    points_new.append((x, y_new))

with open(filename+"_norm", "w") as fout:
    for point in points_new:
        x, y = point
        fout.write(f"{x}\t{y}\n")

