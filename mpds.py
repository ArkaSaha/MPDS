import sys
from os import system

assert len(sys.argv) == 6, "Usage: python mpds.py path-to-graph path-to-output number-of-samples number-of-subgraphs notion-of-density"

graph = sys.argv[1]
output = sys.argv[2]
theta = int(sys.argv[3])
k = int(sys.argv[4])
notion = sys.argv[5]
intake = lambda msg : input(msg) if sys.version_info[0] >= 3 else raw_input(msg)

if notion == "edge":
	executable = "edge/mpds_edge"
	# system("rm -f {ex} && g++ -O3 -std=c++11 -o {ex} {ex}.cpp -lpthread".format(ex = executable))
	system("{} {} {} 40 {} {}".format(executable, graph, theta, k, output))
elif notion == "clique":
	h = intake("Value of h (3/4/5): ")
	executable = "clique/mpds_clique"
	# system("rm -f {ex} && g++ -O3 -std=c++11 -o {ex} {ex}.cpp -lpthread".format(ex = executable))
	system("{} {} {} 40 {} {} {}".format(executable, graph, theta, k, h, output))
elif notion == "pattern":
	method = "heuristic" if intake("Heuristic method (yes /no): ") == "yes" else "approx"
	psi = intake("Pattern (2-star/3-star/c3-star/diamond): ")
	executable = "pattern/mpds_pattern_{}".format(method)
	# system("rm -f {ex} && g++ -O3 -std=c++11 -o {ex} {ex}.cpp -lpthread".format(ex = executable))
	system("{} {} {} 40 {} {} {}".format(executable, graph, theta, k, psi, output))
else:
	assert False, "Invalid notion of density"
