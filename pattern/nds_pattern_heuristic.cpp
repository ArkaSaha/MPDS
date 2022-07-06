#include "util.hpp"

vector< vector<size_t> > core_reduce_diamond(unordered_map< size_t, unordered_set<size_t> >& graph, vector< vector<size_t> >& choose)
{
	unordered_map<size_t, size_t> pat_deg = unordered_map<size_t, size_t>();
	size_t m = 0, n = graph.size(), max_core = 0;
	for (auto& e : graph)
	{
		size_t v = e.first;
		unordered_map<size_t, size_t> paths = unordered_map<size_t, size_t>();
		for (size_t u : e.second)
			for (size_t w : graph[u])
				if (w != v)
				{
					if (paths.find(w) == paths.end())
						paths[w] = 0;
					paths[w]++;
				}
		pat_deg[v] = 0;
		for (auto& p : paths)
			pat_deg[v] += choose[p.second][2];
		m += pat_deg[v];
	}
	m /= 4;
	double max_den = (double) m / n;
	fib_heap heap = fib_heap();
	unordered_map<size_t, handle_t> handles = unordered_map<size_t, handle_t>();
	for (auto& d : pat_deg)
		handles[d.first] = heap.push(node(d.first, d.second));
	vector< pair<size_t, double> > rem = vector< pair<size_t, double> >();
	rem.push_back(make_pair(numeric_limits<int>::max(), max_den));
	while (!heap.empty())
	{
		node nd = heap.top();
		heap.pop();
		size_t v = nd.vertex, d = nd.degree;
		handles.erase(v);
		unordered_map<size_t, size_t> paths = unordered_map<size_t, size_t>(), de = unordered_map<size_t, size_t>();
		for (size_t u : graph[v])
			if (handles.find(u) != handles.end())
				for (size_t w : graph[u])
					if (handles.find(w) != handles.end())
					{
						if (paths.find(w) == paths.end())
							paths[w] = 0;
						paths[w]++;
					}
		for (size_t u : graph[v])
			if (handles.find(u) != handles.end())
			{
				de[u] = 0;
				for (size_t w : graph[u])
					if (handles.find(w) != handles.end())
						de[u] += (paths[w] - 1);
			}
		for (auto& p : paths)
		{
			if (de.find(p.first) == de.end())
				de[p.first] = 0;
			de[p.first] += choose[p.second][2];
			m -= choose[p.second][2];
		}
		for (auto& e : de)
			heap.update(handles[e.first], node(e.first, max((*handles[e.first]).degree - e.second, d)));
		n--;
		double den = (double) m / n;
		rem.push_back(make_pair(v, den));
		if (max_core < d)
		{
			max_core = d;
			max_den = den;
		}
		pat_deg.erase(v);
	}
	vector<size_t> nodes = vector<size_t>();
	vector< vector<size_t> > res = vector< vector<size_t> >();
	while (!rem.empty())
	{
		pair<size_t, double> r = rem.back();
		rem.pop_back();
		if (r.second >= max_den)
			res.push_back(nodes);
		nodes.push_back(r.first);
	}
	return res;
}

vector< vector<size_t> > core_reduce_star(unordered_map< size_t, unordered_set<size_t> >& graph, vector< vector<size_t> >& choose, size_t x)
{
	unordered_map<size_t, size_t> deg = unordered_map<size_t, size_t>(), pat_deg = unordered_map<size_t, size_t>();
	size_t m = 0, n = graph.size(), max_core = 0;
	for (auto& e : graph)
	{
		size_t v = e.first;
		deg[v] = e.second.size();
		pat_deg[v] = choose[deg[v]][x];
		for (size_t u : e.second)
			if (!graph[u].empty())
				pat_deg[v] += choose[graph[u].size() - 1][x - 1];
		m += pat_deg[v];
	}
	m /= (x + 1);
	double max_den = (double) m / n;
	fib_heap heap = fib_heap();
	unordered_map<size_t, handle_t> handles = unordered_map<size_t, handle_t>();
	for (auto& d : pat_deg)
		handles[d.first] = heap.push(node(d.first, d.second));
	vector< pair<size_t, double> > rem = vector< pair<size_t, double> >();
	rem.push_back(make_pair(numeric_limits<int>::max(), max_den));
	while (!heap.empty())
	{
		node nd = heap.top();
		heap.pop();
		size_t v = nd.vertex, d = nd.degree;
		handles.erase(v);
		size_t count = 0;
		for (size_t u : graph[v])
			if (handles.find(u) != handles.end() and deg[u] and deg[v])
			{
				size_t rem = (choose[deg[u] - 1][x - 1] + choose[deg[v] - 1][x - 1]);
				assert((*handles[u]).degree >= d);
				heap.update(handles[u], node(u, max((*handles[u]).degree - rem, d)));
				count += rem;
				for (size_t w : graph[u])
					if (handles.find(w) != handles.end() and w != v and deg[u] >= 2)
					{
						rem = choose[deg[u] - 2][x - 2];
						heap.update(handles[w], node(w, max((*handles[w]).degree - rem, d)));
						count += rem;
					}
			}
		m -= (count / x);
		n--;
		double den = (double) m / n;
		rem.push_back(make_pair(v, den));
		if (max_core < d)
		{
			max_core = d;
			max_den = den;
		}
		for (size_t u : graph[v])
			deg[u]--;
		deg.erase(v);
		pat_deg.erase(v);
	}
	vector<size_t> nodes = vector<size_t>();
	vector< vector<size_t> > res = vector< vector<size_t> >();
	while (!rem.empty())
	{
		pair<size_t, double> r = rem.back();
		rem.pop_back();
		if (r.second >= max_den)
			res.push_back(nodes);
		nodes.push_back(r.first);
	}
	return res;
}

vector< vector<size_t> > core_reduce_cstar(unordered_map< size_t, unordered_set<size_t> >& graph)
{
	unordered_map<size_t, size_t> pat_deg = unordered_map<size_t, size_t>();
	unordered_map< size_t, unordered_set<size_t> > adj = unordered_map< size_t, unordered_set<size_t> >();
	size_t m = 0, max_core = 0;
	vector< set<size_t> > cliques = get_cliques(graph, 3);
	for (size_t i = 0; i < cliques.size(); i++)
	{
		set<size_t> c = cliques[i];
		size_t count = 0;
		for (size_t v : c)
			for (size_t u : graph[v])
				if (c.find(u) == c.end())
				{
					count++;
					if (pat_deg.find(u) == pat_deg.end())
						pat_deg[u] = 0;
					pat_deg[u]++;
				}
		m += count;
		for (size_t v : c)
		{
			if (adj.find(v) == adj.end())
				adj[v] = unordered_set<size_t>();
			adj[v].insert(i);
			if (pat_deg.find(v) == pat_deg.end())
				pat_deg[v] = 0;
			pat_deg[v] += count;
		}
	}
	for (auto& e : graph)
		if (pat_deg.find(e.first) == pat_deg.end())
			pat_deg[e.first] = 0;
	size_t n = pat_deg.size();
	double max_den = (double) m / n;
	fib_heap heap = fib_heap();
	unordered_map<size_t, handle_t> handles = unordered_map<size_t, handle_t>();
	for (auto& d : pat_deg)
		handles[d.first] = heap.push(node(d.first, d.second));
	vector< pair<size_t, double> > rem = vector< pair<size_t, double> >();
	rem.push_back(make_pair(numeric_limits<int>::max(), max_den));
	while (!heap.empty())
	{
		node nd = heap.top();
		heap.pop();
		size_t v = nd.vertex, d = nd.degree;
		handles.erase(v);
		unordered_map<size_t, size_t> changed = unordered_map<size_t, size_t>();
		for (size_t i : adj[v])
		{
			set<size_t> c = cliques[i];
			size_t count = 0;
			for (size_t u : c)
				for (size_t w : graph[u])
					if (handles.find(w) != handles.end() and c.find(w) == c.end())
					{
						count++;
						if (changed.find(w) == changed.end())
							changed[w] = 0;
						changed[w]++;
					}
			for (size_t u : c)
				if (u != v)
				{
					if (changed.find(u) == changed.end())
						changed[u] = 0;
					changed[u] += count;
					adj[u].erase(i);
				}
			m -= count;
		}
		for (size_t u : graph[v])
			if (handles.find(u) != handles.end())
				for (size_t i : adj[u])
				{
					for (size_t w : cliques[i])
					{
						if (changed.find(w) == changed.end())
							changed[w] = 0;
						changed[w]++;
					}
					m--;
				}
		for (auto& c : changed)
			heap.update(handles[c.first], node(c.first, max((*handles[c.first]).degree - c.second, d)));
		n--;
		double den = (double) m / n;
		rem.push_back(make_pair(v, den));
		if (max_core < d)
		{
			max_core = d;
			max_den = den;
		}
		pat_deg.erase(v);
		adj.erase(v);
	}
	vector<size_t> nodes = vector<size_t>();
	vector< vector<size_t> > res = vector< vector<size_t> >();
	while (!rem.empty())
	{
		pair<size_t, double> r = rem.back();
		rem.pop_back();
		if (r.second >= max_den)
			res.push_back(nodes);
		nodes.push_back(r.first);
	}
	return res;
}

vector< vector<size_t> > pattern_exact(unordered_map< size_t, unordered_set<size_t> >& graph, string pattern, vector< vector<size_t> >& choose)
{
	if (pattern == "c3-star")
		return core_reduce_cstar(graph);
	else if (pattern == "diamond")
		return core_reduce_diamond(graph, choose);
	else
		return core_reduce_star(graph, choose, pattern[0] - '0');
}

struct thread_data
{
	int samples;
	string pattern;
	vector< vector<size_t> > choose;
	vector< tuple<size_t, size_t, double> > edges;
	vector< unordered_set<size_t> > candidates;
};

void* helper(void* arg)
{
	thread_data* tdata = (thread_data*) arg;
	tdata->candidates = vector< unordered_set<size_t> >();
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
		vector< vector<size_t> > cand = pattern_exact(graph, tdata->pattern, tdata->choose);
		unordered_set<size_t> nodes = unordered_set<size_t>();
		for (vector<size_t> vs : cand)
			for (size_t v : vs)
				nodes.insert(v);
		tdata->candidates.push_back(nodes);
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
		for (unordered_set<size_t> c : tdata[i].candidates)
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
