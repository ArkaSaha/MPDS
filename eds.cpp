#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <boost/heap/fibonacci_heap.hpp>
#define BOOST_DISABLE_ASSERTS
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/graph_utility.hpp>

using namespace std;
using namespace boost::heap;
using Traits = boost::adjacency_list_traits < boost::vecS, boost::vecS, boost::directedS >;
using Vertex = Traits::vertex_descriptor;
using Edge = Traits::edge_descriptor;
using Graph = boost::adjacency_list < boost::vecS, boost::vecS, boost::directedS, boost::property < boost::vertex_index_t, size_t >, boost::property < boost::edge_capacity_t, double, boost::property < boost::edge_residual_capacity_t, double, boost::property < boost::edge_reverse_t, Edge > > > >;

void core_reduce(unordered_map< size_t, unordered_map<size_t, double> >& adj)
{
	unordered_map<size_t, double> deg = unordered_map<size_t, double>(), core = unordered_map<size_t, double>();
	size_t n = adj.size();
	double m = 0;
	for (auto& e : adj)
		for (auto& f : e.second)
		{
			m += f.second;
			deg[e.first] += f.second;
		}
	m /= 2;
	double max_den = m / n;
	struct node
	{
		size_t vertex;
		double degree;
		node(const size_t& v, double d) : vertex(v), degree(d) {}
	};
	struct compare_node
	{
		bool operator()(const node& n1, const node& n2) const
		{
			return n1.degree > n2.degree;
		}
	};
	using handle_t = fibonacci_heap< node, compare<compare_node> >::handle_type;
	fibonacci_heap< node, compare<compare_node> > heap = fibonacci_heap< node, compare<compare_node> >();
	unordered_map<size_t, handle_t> handles = unordered_map<size_t, handle_t>();
	for (auto& d : deg)
		handles[d.first] = heap.push(node(d.first, d.second));
	while (!heap.empty())
	{
		node nd = heap.top();
		heap.pop();
		size_t v = nd.vertex;
		double d = nd.degree;
		handles.erase(v);
		core[v] = d;
		for (auto& e : adj[v])
			if (handles.find(e.first) != handles.end())
			{
				heap.update(handles[e.first], node(e.first, max((*handles[e.first]).degree - e.second, d)));
				m -= e.second;
			}
		n--;
		double den = m / n;
		if (n and den > max_den)
			max_den = den;
		deg.erase(v);
	}
	for (auto& d : core)
	{
		size_t v = d.first;
		if (d.second >= max_den)
		{
			unordered_set<size_t> rem = unordered_set<size_t>();
			for (auto& e : adj[v])
				if (core[e.first] < max_den)
					rem.insert(e.first);
			for (size_t u : rem)
				adj[v].erase(u);
		}
		else
			adj.erase(v);
	}
}

unordered_map<size_t, double> kclistpp(unordered_map< size_t, unordered_map<size_t, double> >& adj, size_t T, unordered_map< size_t, unordered_map<size_t, double> >& alpha)
{
	unordered_map<size_t, double> r = unordered_map<size_t, double>();
	for (auto& a : alpha)
	{
		r[a.first] = 0;
		for (auto& b : a.second)
			r[a.first] += b.second;
	}
	for (size_t t = 1; t <= T; t++)
	{
		double gamma = 1.0 / (t + 1);
		for (auto& e : alpha)
			for (auto& f : e.second)
				alpha[e.first][f.first] *= (1 - gamma);
		for (auto& e : r)
			r[e.first] *= (1 - gamma);
		for (auto& e : adj)
		{
			size_t u = e.first;
			for (auto& f : e.second)
			{
				size_t v = f.first;
				double c = f.second;
				if (u < v)
				{
					if (r[v] < r[u])
					{
						alpha[v][u] += (gamma * c);
						r[v] += (gamma * c);
					}
					else
					{
						alpha[u][v] += (gamma * c);
						r[u] += (gamma * c);
					}
				}
			}
		}
	}
	return r;
}

pair<size_t, double> extract_densest(unordered_map< size_t, unordered_map<size_t, double> >& adj, vector<size_t>& nodes, unordered_map< size_t, unordered_map<size_t, double> >& alpha, unordered_map<size_t, double>& r)
{
	vector<size_t> nag = vector<size_t>();
	vector<double> val = vector<double>();
	vector<double> tentative = vector<double>(nodes.size(), 0);
	unordered_map<size_t, size_t> level = unordered_map<size_t, size_t>();
	vector< pair<size_t, double> > tmp = vector< pair<size_t, double> >();
	for (size_t v : nodes)
		tmp.push_back(make_pair(v, r[v]));
	struct { bool operator()(pair<size_t, double> p, pair<size_t, double> q) {return p.second > q.second;} } comp;
	sort(tmp.begin(), tmp.end(), comp);
	for (size_t i = 0; i < nodes.size(); i++)
		nodes[i] = tmp[i].first;
	for (size_t i = 0; i < nodes.size(); i++)
		level[nodes[i]] = i;
	for (auto& e : adj)
	{
		size_t u = e.first;
		for (auto& f : e.second)
		{
			size_t v = f.first;
			double c = f.second;
			if (u < v)
			{
				if (level[v] > level[u])
					tentative[level[v]] += c;
				else
					tentative[level[u]] += c;
			}
		}
	}
	size_t j = 0;
	val.push_back(tentative[0]);
	nag.push_back(1);
	for (size_t i = 1; i < nodes.size(); i++)
	{
		j++;
		val.push_back(tentative[i]);
		nag.push_back(1);
		while (j > 0 and val[j] >= val[j - 1] * (1 - 1e-6))
		{
			val[j - 1] = (nag[j] * val[j] + nag[j - 1] * val[j - 1]) / (nag[j] + nag[j - 1]);
			nag[j - 1] += nag[j];
			val.pop_back();
			nag.pop_back();
			j--;
		}
	}
	for (size_t k = 0, i = 0; k < j + 1; ++k)
		for (size_t l = 0; l < nag[k]; ++l, ++i)
			level[nodes[i]] = k;
	for (auto& e : adj)
	{
		size_t u = e.first;
		for (auto& f : e.second)
		{
			size_t v = f.first;
			if (u < v)
			{
				size_t max_level = max(level[u], level[v]), max_level_cnt = (level[u] == level[v]) ? 2 : 1;
				double sum = 0;
				if (level[u] < max_level)
				{
					sum += alpha[u][v];
					r[u] -= alpha[u][v];
					alpha[u][v] = 0;
				}
				if (level[v] < max_level)
				{
					sum += alpha[v][u];
					r[v] -= alpha[v][u];
					alpha[v][u] = 0;
				}
				if (level[u] == max_level)
				{
					r[u] += (sum / max_level_cnt);
					alpha[u][v] += (sum / max_level_cnt);
				}
				if (level[v] == max_level)
				{
					r[v] += (sum / max_level_cnt);
					alpha[v][u] += (sum / max_level_cnt);
				}
			}
		}
	}
	return make_pair(nag[0], val[0]);
}

void insert_edge(Graph& g, Vertex& v1, Vertex& v2, double c1, double c2 = 0)
{
	Edge e1 = boost::add_edge(v1, v2, g).first;
	Edge e2 = boost::add_edge(v2, v1, g).first;
	put(boost::edge_capacity, g, e1, c1);
	put(boost::edge_capacity, g, e2, c2);
	get(boost::edge_reverse, g)[e1] = e2;
	get(boost::edge_reverse, g)[e2] = e1;
}

bool densest_maxflow(unordered_map< size_t, unordered_map<size_t, double> >& adj, vector<size_t>& nodes, size_t n, double den)
{
    Graph g = Graph();
    vector<Vertex> R = vector<Vertex>();
	Vertex s = boost::add_vertex(g), t = boost::add_vertex(g);
	unordered_map<size_t, size_t> ids = unordered_map<size_t, size_t>();
	for (size_t i = 0; i < nodes.size(); i++)
		ids[nodes[i]] = i;
	for (size_t i = 0; i < n; i++)
		R.push_back(boost::add_vertex(g));
	for (auto& e : adj)
	{
		size_t u = e.first;
		if (ids[u] < n)
		{
			double deg = 0;
			for (auto& f : e.second)
			{
				size_t v = f.first;
				double c = f.second;
				if (ids[v] < n)
				{
					deg += c;
					if (u < v)
						insert_edge(g, R[ids[u]], R[ids[v]], c, c);
				}
			}
			insert_edge(g, s, R[ids[u]], deg);
			insert_edge(g, R[ids[u]], t, 2 * den);
		}
	}
	return boost::push_relabel_max_flow(g, s, t) >= 2 * n * den * (1 - 1e-6);
}

pair<size_t, double> maximum_density(unordered_map< size_t, unordered_map<size_t, double> >& adj, vector<size_t>& nodes)
{
	unordered_map< size_t, unordered_map<size_t, double> > alpha = unordered_map< size_t, unordered_map<size_t, double> >();
	for (auto& e : adj)
	{
		size_t u = e.first;
		if (alpha.find(u) == alpha.end())
			alpha[u] = unordered_map<size_t, double>();
		for (auto& f : e.second)
			alpha[u][f.first] = f.second / 2;
	}
	for (size_t t = 1; ; t *= 2)
	{
		unordered_map<size_t, double> r = kclistpp(adj, t, alpha);
		pair<size_t, double> res = extract_densest(adj, nodes, alpha, r);
		double prefix_min_rho = r[nodes[0]], suffix_max_rho = -1;
		for (size_t i = 1; i < res.first; ++i)
			if (prefix_min_rho > r[nodes[i]])
				prefix_min_rho = r[nodes[i]];
		for (size_t i = nodes.size() - 1; i >= res.first; --i)
			if (suffix_max_rho < r[nodes[i]])
				suffix_max_rho = r[nodes[i]];
		if (prefix_min_rho * (1 - 1e-6) > suffix_max_rho and densest_maxflow(adj, nodes, res.first, res.second))
			return res;
	}
	return make_pair(0, 0);
}

void eds(char* output, unordered_map< size_t, unordered_map<size_t, double> >& adj)
{
	vector<size_t> nodes = vector<size_t>();
	timespec begin, end;
	clock_gettime(CLOCK_MONOTONIC, &begin);
	// core_reduce(adj);
	for (auto& e : adj)
		nodes.push_back(e.first);
	pair<size_t, double> res = maximum_density(adj, nodes);
	clock_gettime(CLOCK_MONOTONIC, &end);
	cout << "Running time : " << (end.tv_sec - begin.tv_sec) + (end.tv_nsec - begin.tv_nsec) / pow(10, 9) << " seconds" << endl;
	ofstream fout(output);
	for (size_t i = 0; i < res.first; i++)
		fout << nodes[i] << " ";
	fout << res.second << endl;
	fout.close();
}

int main(int argc, char* argv[])
{
	srand(time(NULL));
	if (argc != 3)
	{
		cerr << "Usage : ./edge path-to-graph path-to-output" << endl;
		return EXIT_FAILURE;
	}
	size_t n, u, v;
	double p;
	unordered_map< size_t, unordered_map<size_t, double> > adj = unordered_map< size_t, unordered_map<size_t, double> >();
	ifstream fin(argv[1]);
	fin >> n;
	while (fin >> u >> v >> p)
		if (u != v)
		{
			if (adj.find(u) == adj.end())
				adj[u] = unordered_map<size_t, double>();
			adj[u][v] = p;
			if (adj.find(v) == adj.end())
				adj[v] = unordered_map<size_t, double>();
			adj[v][u] = p;
		}
	fin.close();
	eds(argv[2], adj);
	return EXIT_SUCCESS;
}
