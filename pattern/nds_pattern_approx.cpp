#include "pattern.hpp"

vector<size_t> pattern_exact(unordered_map< size_t, unordered_set<size_t> >& graph, string pattern, vector< vector<size_t> >& choose)
{
	if (pattern == "c3-star")
		core_reduce_cstar(graph);
	else if (pattern == "diamond")
		core_reduce_diamond(graph, choose);
	else
		core_reduce_star(graph, choose, pattern[0] - '0');
	unordered_map<string, size_t> instances = get_patterns(graph, pattern);
	vector<size_t> nodes = vector<size_t>();
	for (auto& e : graph)
		nodes.push_back(e.first);
	size_t h = pattern == "2-star" ? 3 : 4, factor = pattern == "c3-star" ? 12 : (pattern == "diamond" ? 3 : pattern[0] - '0' + 1);
	pair<size_t, double> res = maximum_density(instances, h, nodes, factor);
	nodes.resize(res.first);
	return nodes;
}

struct thread_data
{
	int samples;
	string pattern;
	vector< vector<size_t> > choose;
	vector< tuple<size_t, size_t, double> > edges;
	vector< vector<size_t> > candidates;
};

void* helper(void* arg)
{
	thread_data* tdata = (thread_data*) arg;
	tdata->candidates = vector< vector<size_t> >();
	for (int i = 1; i <= tdata->samples; i++)
	{
		unordered_map< size_t, unordered_set<size_t> > graph = unordered_map< size_t, unordered_set<size_t> >();
		for (auto e : tdata->edges)
		{
			double r = (double) rand() / RAND_MAX;
			if (r < get<2>(e))
			{
				size_t u = get<0>(e), v = get<1>(e);
				if (graph.find(u) == graph.end())
					graph[u] = unordered_set<size_t>();
				graph[u].insert(v);
				if (graph.find(v) == graph.end())
					graph[v] = unordered_set<size_t>();
				graph[v].insert(u);
			}
		}
		tdata->candidates.push_back(pattern_exact(graph, tdata->pattern, tdata->choose));
	}
	pthread_exit(NULL);
}

void mpds(int theta, int num_threads, string pattern, char* output, vector< tuple<size_t, size_t, double> > edges)
{
	int per_thread = ceil((double) theta / num_threads);
	timespec begin, end;
	vector<thread_data> tdata = vector<thread_data>(num_threads);
	vector<pthread_t> tid = vector<pthread_t>(num_threads);
	unordered_map<size_t, size_t> deg = unordered_map<size_t, size_t>();
	for (tuple<size_t, size_t, double> e : edges)
	{
		size_t u = get<0>(e), v = get<1>(e);
		if (deg.find(u) == deg.end())
			deg[u] = 0;
		deg[u]++;
		if (deg.find(v) == deg.end())
			deg[v] = 0;
		deg[v]++;
	}
	size_t max_deg = 0;
	for (auto d : deg)
		if (d.second > max_deg)
			max_deg = d.second;
	vector< vector<size_t> > choose = vector< vector<size_t> >(max_deg + 1, vector<size_t>(4, 0));
	for (size_t i = 0; i <= max_deg; i++)
		choose[i][0] = 1;
	for (size_t j = 1; j < 4; j++)
	{
		choose[j][j] = 1;
		for (size_t i = j; i < max_deg; i++)
			choose[i + 1][j] = choose[i][j] + choose[i][j - 1];
	}
	clock_gettime(CLOCK_MONOTONIC, &begin);
	for (int i = 0; i < num_threads; i++)
	{
		tdata[i].edges = edges;
		tdata[i].samples = per_thread;
		tdata[i].pattern = pattern;
		tdata[i].choose = choose;
		if (pthread_create(&tid[i], NULL, helper, (void*)(&tdata[i])))
			cerr << "ERROR : pthread_create" << endl;
	}
	for (int i = 0; i < num_threads; i++)
		pthread_join(tid[i], NULL);
	clock_gettime(CLOCK_MONOTONIC, &end);
	cout << "Running time : " << (end.tv_sec - begin.tv_sec) + (end.tv_nsec - begin.tv_nsec) / pow(10, 9) << " seconds" << endl;
	ofstream fout(output);
	for (int i = 0; i < num_threads; i++)
	{
		for (vector<size_t> c : tdata[i].candidates)
		{
			for (size_t v : c)
				fout << v << " ";
			fout << endl;
		}
	}
	fout << endl;
	fout.close();
}

int main(int argc, char* argv[])
{
	srand(time(NULL));
	if (argc != 6)
	{
		cerr << "Usage : ./edge path-to-graph number-of-samples number-of-threads pattern path-to-output" << endl;
		return EXIT_FAILURE;
	}
	size_t n, u, v;
	double p;
	vector< tuple<size_t, size_t, double> > edges = vector< tuple<size_t, size_t, double> >();
	ifstream fin(argv[1]);
	fin >> n;
	while (fin >> u >> v >> p)
		edges.push_back(make_tuple(u, v, p));
	fin.close();
	mpds(atoi(argv[2]), atoi(argv[3]), argv[4], argv[5], edges);
	return EXIT_SUCCESS;
}
