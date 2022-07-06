import sys
from os import system, popen

assert len(sys.argv) == 6, "Usage: python nds.py path-to-graph path-to-output number-of-samples containment-probability-threshold notion-of-density"

graph = sys.argv[1]
output = sys.argv[2]
theta = int(sys.argv[3])
sigma = float(sys.argv[4])
notion = sys.argv[5]
intake = lambda msg : input(msg) if sys.version_info[0] >= 3 else raw_input(msg)

runtime = 0
if notion == "edge":
	executable = "edge/nds_edge"
	system("rm -f {ex} && g++ -O3 -lpthread -std=c++11 -o {ex} {ex}.cpp".format(ex = executable))
	runtime = float(popen("{} {} {} 40 input.log".format(executable, graph, theta)).read()[15:-9])
elif notion == "clique":
	h = intake("Value of h (3/4/5): ")
	executable = "clique/nds_clique"
	system("rm -f {ex} && g++ -O3 -lpthread -std=c++11 -o {ex} {ex}.cpp".format(ex = executable))
	runtime = float(popen("{} {} {} 40 {} input.log".format(executable, graph, theta, h)).read()[15:-9])
elif notion == "pattern":
	method = "heuristic" if intake("Heuristic method (yes /no): ") == "yes" else "approx"
	psi = intake("Pattern (2-star/3-star/c3-star/diamond): ")
	executable = "pattern/nds_pattern_{}".format(method)
	system("rm -f {ex} && g++ -O3 -lpthread -std=c++11 -o {ex} {ex}.cpp".format(ex = executable))
	runtime = float(popen("{} {} {} 40 {} input.log".format(executable, graph, theta, psi)).read()[15:-9])
else:
	assert False, "Invalid notion of density"
start = time()
system("./fpgrowth -m0tms{} input.log {}".format(sigma * (theta - 0.5) / (theta + 1), output))
end = time()
print("Running time : ({} + {}) seconds = {} seconds".format(runtime, end - start, runtime + end - start))
