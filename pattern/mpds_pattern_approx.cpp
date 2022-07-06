#include "pattern.hpp"

vector< vector<size_t> > pattern_exact(unordered_map< size_t, unordered_set<size_t> >& graph, size_t n, string pattern, vector< vector<size_t> >& choose)
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
	unordered_map< size_t, vector<size_t> > graph_new = construct(n, h, instances, res.second);
	unordered_map< size_t, unordered_set<size_t> > scc_nodes = unordered_map< size_t, unordered_set<size_t> >();
	unordered_map< size_t, unordered_set<size_t> > desc = unordered_map< size_t, unordered_set<size_t> >();
	unordered_map< size_t, unordered_set<size_t> > anc = unordered_map< size_t, unordered_set<size_t> >();
	unordered_map<size_t, unordered_set<size_t> > closed_nodes = unordered_map<size_t, unordered_set<size_t> >();
	compute_scc(graph_new, n, scc_nodes, desc, closed_nodes);
	unordered_set<size_t> c = unordered_set<size_t>();
	for (auto e : desc)
	{
		size_t u = e.first;
		c.insert(u);
		if (anc.find(u) == anc.end())
			anc[u] = unordered_set<size_t>();
		for (size_t v : e.second)
		{
			if (anc.find(v) == anc.end())
				anc[v] = unordered_set<size_t>();
			anc[v].insert(u);
		}
	}
	vector< vector<size_t> > subgraphs = vector< vector<size_t> >();
	enumerate_densest(unordered_set<size_t>(), c, numeric_limits<int>::max(), n, desc, anc, scc_nodes, closed_nodes, subgraphs);
	nodes.resize(res.first);
	subgraphs.push_back(nodes);
	return subgraphs;
}

struct thread_data
{
	int samples;
	size_t n;
	string pattern;
	vector< vector<size_t> > choose;
	vector< tuple<size_t, size_t, double> > edges;
	unordered_map<string,size_t> freq;
};

void* helper(void* arg)
{
	thread_data* tdata = (thread_data*) arg;
	tdata->freq = unordered_map<string,size_t>();
	for (int i = 1; i <= tdata->samples; i++)
	{
		unordered_map< size_t, unordered_set<size_t> > graph = unordered_map< size_t, unordered_set<size_t> >();
		for (auto e : tdata->edges)
			if ((double) rand() / RAND_MAX < get<2>(e))
			{
				size_t u = get<0>(e), v = get<1>(e);
				if (graph.find(u) == graph.end())
					graph[u] = unordered_set<size_t>();
				graph[u].insert(v);
				if (graph.find(v) == graph.end())
					graph[v] = unordered_set<size_t>();
				graph[v].insert(u);
			}
		vector< vector<size_t> > cand = pattern_exact(graph, tdata->n, tdata->pattern, tdata->choose);
		unordered_set<string> c = unordered_set<string>();
		for (vector<size_t> v : cand)
		{
			string h = code(v);
			if (c.find(h) == c.end())
			{
				if (tdata->freq.find(h) == tdata->freq.end())
					tdata->freq[h] = 0;
				tdata->freq[h]++;
				c.insert(h);
			}
		}
	}
	pthread_exit(NULL);
}

void mpds(int theta, int num_threads, size_t k, string pattern, char* output, size_t n, vector< tuple<size_t, size_t, double> > edges)
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
		tdata[i].n = n;
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
	unordered_map<string,size_t> freq = unordered_map<string,size_t>();
	for (int i = 0; i < num_threads; i++)
	{
		for (auto f : tdata[i].freq)
		{
			if (freq.find(f.first) == freq.end())
				freq[f.first] = 0;
			freq[f.first] += f.second;
		}
	}
	vector< pair<vector<size_t>,size_t> > cv = vector< pair<vector<size_t>,size_t> >();
	for (auto f : freq)
		cv.push_back(make_pair(decode(f.first), f.second));
	struct { bool operator() (pair<vector<size_t>,size_t> x, pair<vector<size_t>,size_t> y) { return x.second > y.second; } } comp;
	sort(cv.begin(), cv.end(), comp);
	k = min(k, cv.size());
	ofstream fout(output);
	for (size_t i = 0; i < k; i++)
	{
		for (size_t v : cv[i].first)
			fout << v << " ";
		double prob = (double) cv[i].second / theta;
		fout << "\t" << prob << endl;
	}
	fout.close();
}

int main(int argc, char* argv[])
{
	srand(time(NULL));
	if (argc != 7)
	{
		cerr << "Usage : ./edge path-to-graph number-of-samples number-of-threads number-of-subgraphs pattern path-to-output" << endl;
		return EXIT_FAILURE;
	}
	size_t n, u, v;
	double p;
	vector< tuple<size_t, size_t, double> > edges = vector< tuple<size_t, size_t, double> >();
	ifstream fin(argv[1]);
	fin >> n;
	while (fin >> u >> v >> p)
		if (u != v)
			edges.push_back(make_tuple(u, v, p));
	fin.close();
	mpds(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), argv[5], argv[6], n, edges);
	return EXIT_SUCCESS;
}
