#include "edge.hpp"

vector< vector<size_t> > edge_exact(unordered_map< size_t, unordered_set<size_t> >& adj, size_t n)
{
	double max_den = 0;
	core_reduce(adj, max_den);
	vector<size_t> nodes = vector<size_t>();
	for (auto& e : adj)
		nodes.push_back(e.first);
	pair<size_t, double> res = maximum_density(adj, nodes);
	unordered_map< size_t, vector<size_t> > adj_residual = construct(n, adj, res.second);
	unordered_map< size_t, unordered_set<size_t> > scc_nodes = unordered_map< size_t, unordered_set<size_t> >();
	unordered_map< size_t, unordered_set<size_t> > desc = unordered_map< size_t, unordered_set<size_t> >();
	unordered_map< size_t, unordered_set<size_t> > anc = unordered_map< size_t, unordered_set<size_t> >();
	unordered_map<size_t, unordered_set<size_t> > closed_nodes = unordered_map<size_t, unordered_set<size_t> >();
	compute_scc(adj_residual, n, scc_nodes, desc, closed_nodes);
	unordered_set<size_t> c = unordered_set<size_t>();
	for (auto& e : desc)
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
	enumerate_densest(unordered_set<size_t>(), c, INT_MAX, n, desc, anc, scc_nodes, closed_nodes, subgraphs);
	nodes.resize(res.first);
	subgraphs.push_back(nodes);
	return subgraphs;
}

void helper(vector< tuple<size_t, size_t, double> >& edges, size_t n, size_t pos, vector<size_t>& selected, double prob, unordered_map<string, double>& freq)
{
	if (pos == edges.size())
	{
		if (selected.empty())
			return;
		unordered_map< size_t, unordered_set<size_t> > adj = unordered_map< size_t, unordered_set<size_t> >();
		for (size_t i : selected)
		{
			tuple<size_t, size_t, double> e = edges[i];
			size_t u = get<0>(e), v = get<1>(e);
			if (adj.find(u) == adj.end())
				adj[u] = unordered_set<size_t>();
			adj[u].insert(v);
			if (adj.find(v) == adj.end())
				adj[v] = unordered_set<size_t>();
			adj[v].insert(u);
		}
		vector< vector<size_t> > cand = edge_exact(adj, n);
		unordered_set<string> c = unordered_set<string>();
		for (vector<size_t>& v : cand)
		{
			string h = code(v);
			if (c.find(h) == c.end())
			{
				if (freq.find(h) == freq.end())
					freq[h] = 0;
				freq[h] += prob;
				c.insert(h);
			}
		}
	}
	else
	{
		double pr = get<2>(edges[pos]);
		helper(edges, n, pos + 1, selected, prob * (1 - pr), freq);
		selected.push_back(pos);
		helper(edges, n, pos + 1, selected, prob * pr, freq);
		selected.pop_back();
	}
}

void mpds(size_t k, char* output, size_t n, vector< tuple<size_t, size_t, double> >& edges)
{
	timespec begin, end;
	unordered_map<string, double> freq = unordered_map<string, double>();
	vector<size_t> selected = vector<size_t>();
	clock_gettime(CLOCK_MONOTONIC, &begin);
	helper(edges, n, 0, selected, 1, freq);
	clock_gettime(CLOCK_MONOTONIC, &end);
	cout << "Running time : " << (end.tv_sec - begin.tv_sec) + (end.tv_nsec - begin.tv_nsec) / pow(10, 9) << " seconds" << endl;
	vector< pair<vector<size_t>, double> > cv = vector< pair<vector<size_t>, double> >();
	for (auto& f : freq)
		cv.push_back(make_pair(decode(f.first), f.second));
	struct { bool operator() (pair<vector<size_t>, double> x, pair<vector<size_t>, double> y) { return x.second > y.second; } } comp;
	sort(cv.begin(), cv.end(), comp);
	k = min(k, cv.size());
	ofstream fout(output);
	for (size_t i = 0; i < k; i++)
	{
		for (size_t v : cv[i].first)
			fout << v << " ";
		fout << "\t" << cv[i].second << endl;
	}
	fout.close();
}

int main(int argc, char* argv[])
{
	srand(time(NULL));
	if (argc != 4)
	{
		cerr << "Usage : ./edge path-to-graph number-of-subgraphs path-to-output" << endl;
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
	mpds(atoi(argv[2]), argv[3], n + 1, edges);
	return EXIT_SUCCESS;
}
