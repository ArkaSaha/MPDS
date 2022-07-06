#include "edge.hpp"
#include "memoryusage.hpp"

vector<size_t> edge_exact(unordered_map< size_t, unordered_set<size_t> >& adj)
{
	double max_den = 0;
	core_reduce(adj, max_den);
	vector<size_t> nodes = vector<size_t>();
	for (auto& e : adj)
		nodes.push_back(e.first);
	pair<size_t, double> res = maximum_density(adj, nodes);
	nodes.resize(res.first);
	return nodes;
}

struct thread_data
{
	int samples;
	vector< tuple<size_t, size_t, double> > edges;
	vector< vector<size_t> > candidates;
};

void* helper(void* arg)
{
	thread_data* tdata = (thread_data*) arg;
	tdata->candidates = vector< vector<size_t> >();
	for (int i = 1; i <= tdata->samples; i++)
	{
		unordered_map< size_t, unordered_set<size_t> > adj = unordered_map< size_t, unordered_set<size_t> >();
		for (auto& e : tdata->edges)
			if ((double) rand() / RAND_MAX < get<2>(e))
			{
				size_t u = get<0>(e), v = get<1>(e);
				if (u != v)
				{
					if (adj.find(u) == adj.end())
						adj[u] = unordered_set<size_t>();
					adj[u].insert(v);
					if (adj.find(v) == adj.end())
						adj[v] = unordered_set<size_t>();
					adj[v].insert(u);
				}
			}
		tdata->candidates.push_back(edge_exact(adj));
	}
	pthread_exit(NULL);
}

void mpds(int theta, int num_threads, char* output, vector< tuple<size_t, size_t, double> >& edges)
{
	int per_thread = ceil((double) theta / num_threads);
	timespec begin, end;
	vector<thread_data> tdata = vector<thread_data>(num_threads);
	vector<pthread_t> tid = vector<pthread_t>(num_threads);
	clock_gettime(CLOCK_MONOTONIC, &begin);
	for (int i = 0; i < num_threads; i++)
	{
		tdata[i].edges = edges;
		tdata[i].samples = per_thread;
		if (pthread_create(&tid[i], NULL, helper, (void*)(&tdata[i])))
			cerr << "ERROR : pthread_create" << endl;
	}
	for (int i = 0; i < num_threads; i++)
		pthread_join(tid[i], NULL);
	clock_gettime(CLOCK_MONOTONIC, &end);
	cout << "Running time : " << (end.tv_sec - begin.tv_sec) + (end.tv_nsec - begin.tv_nsec) / pow(10, 9) << " seconds" << endl;
	disp_mem_usage("");
	ofstream fout(output);
	for (int i = 0; i < num_threads; i++)
		for (vector<size_t>& c : tdata[i].candidates)
		{
			for (size_t v : c)
				fout << v << " ";
			fout << endl;
		}
	fout << endl;
	fout.close();
}

int main(int argc, char* argv[])
{
	srand(time(NULL));
	if (argc != 5)
	{
		cerr << "Usage : ./nucleus_edge path-to-graph number-of-samples number-of-threads path-to-output" << endl;
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
	mpds(atoi(argv[2]), atoi(argv[3]), argv[4], edges);
	return EXIT_SUCCESS;
}
