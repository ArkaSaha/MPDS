#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <climits>
#include <ctime>
#include <algorithm>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <pthread.h>
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

struct node
{
	size_t vertex, degree;
	node(const size_t& v, size_t d) : vertex(v), degree(d) {}
};
struct compare_node
{
	bool operator()(const node& n1, const node& n2) const
	{
		return n1.degree > n2.degree;
	}
};
using fib_heap = fibonacci_heap< node, compare<compare_node> >;
using handle_t = fib_heap::handle_type;

void listing(size_t l, vector<size_t>& nodes, set<size_t> c, vector< set<size_t> >& cliques, unordered_map< size_t, vector<size_t> >& dag, unordered_map<size_t, size_t>& deg, unordered_map<size_t, size_t>& label)
{
	if (l == 2)
		for (size_t u : nodes)
			for (size_t i = 0; i < deg[u]; i++)
			{
				size_t v = dag[u][i];
				set<size_t> cc = c;
				cc.insert(u);
				cc.insert(v);
				cliques.push_back(cc);
			}
	else
		for (size_t u : nodes)
		{
			vector<size_t> nodes_new = vector<size_t>();
			for (size_t v : dag[u])
				if (label[v] == l)
				{
					label[v] = l - 1;
					nodes_new.push_back(v);
				}
			for (size_t v : nodes_new)
			{
				size_t index = 0;
				if (!dag[v].empty())
				{
					for (size_t i = dag[v].size() - 1; i > index; i--)
						if (label[dag[v][i]] == l - 1)
						{
							while (index < i and label[dag[v][index]] == l - 1)
								index++;
							if (label[dag[v][index]] != l - 1)
								swap(dag[v][i], dag[v][index]);
						}
					if (label[dag[v][index]] == l - 1)
						index++;
				}
				deg[v] = index;
			}
			c.insert(u);
			listing(l - 1, nodes_new, c, cliques, dag, deg, label);
			c.erase(u);
			for (size_t v : nodes_new)
				label[v] = l;
		}
}

vector< set<size_t> > get_cliques(unordered_map< size_t, unordered_set<size_t> >& adj, size_t h)
{
	unordered_map<size_t, size_t> vert = unordered_map<size_t, size_t>(), core_rev = unordered_map<size_t, size_t>(), deg = unordered_map<size_t, size_t>(), pos = unordered_map<size_t, size_t>();
	unordered_map<size_t, size_t> order = unordered_map<size_t, size_t>();
	size_t md = 0;
	for (auto& e : adj)
	{
		deg[e.first] = e.second.size();
		if (md < deg[e.first])
			md = deg[e.first];
	}
	vector<size_t> bin = vector<size_t>(md + 1, 0);
	for (auto& e : adj)
		bin[deg[e.first]]++;
	size_t start = 1;
	for (size_t d = 0; d <= md; d++)
	{
		size_t num = bin[d];
		bin[d] = start;
		start += num;
	}
	for (auto& d : deg)
	{
		pos[d.first] = bin[d.second];
		vert[pos[d.first]] = d.first;
		bin[d.second]++;
	}
	for (size_t d = md; d >= 1; d--)
		bin[d] = bin[d - 1];
	bin[0] = 1;
	for (size_t i = 1; i <= adj.size(); i++)
	{
		size_t v = vert[i];
		for (size_t u : adj[v])
			if (deg[u] > deg[v])
			{
				size_t du = deg[u], pu = pos[u], pw = bin[du], w = vert[pw];
				if (u != w)
				{
					pos[u] = pw;	vert[pu] = w;
					pos[w] = pu;	vert[pw] = u;
				}
				bin[du]++;
				deg[u]--;
			}
		core_rev[adj.size() - i] = v;
	}
	for (size_t i = 0; i < adj.size(); i++)
		order[core_rev[i]] = i + 1;
	unordered_map< size_t, vector<size_t> > dag = unordered_map< size_t, vector<size_t> >();
	deg.clear();
	vector< pair<size_t, size_t> > tmp = vector< pair<size_t, size_t> >();
	unordered_map<size_t, size_t> label = unordered_map<size_t, size_t>();
	for (auto& e : adj)
	{
		size_t u = e.first;
		if (dag.find(u) == dag.end())
			dag[u] = vector<size_t>();
		if (deg.find(u) == deg.end())
			deg[u] = 0;
		for (size_t v : e.second)
			if (order[u] < order[v])
			{
				dag[u].push_back(v);
				deg[u]++;
			}
		tmp.push_back(make_pair(u, order[u]));
		label[u] = h;
	}
	struct { bool operator()(pair<size_t, size_t> p, pair<size_t, size_t> q) {return p.second < q.second;} } comp;
	sort(tmp.begin(), tmp.end(), comp);
	vector<size_t> nodes = vector<size_t>();
	while (!tmp.empty())
	{
		nodes.push_back(tmp.back().first);
		tmp.pop_back();
	}
	vector< set<size_t> > cliques = vector< set<size_t> >();
	listing(h, nodes, set<size_t>(), cliques, dag, deg, label);
	return cliques;
}

void pava(vector<size_t>& nodes, vector<size_t>& nag, vector<double>& val, vector<size_t>& tentative, unordered_map<size_t, size_t>& level)
{
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
}

bool densest_goldberg(vector<size_t>& nodes, size_t n, size_t h, unordered_map<size_t, double>& r, double den, size_t factor)
{
	vector< pair<size_t, double> > tmp = vector< pair<size_t, double> >();
	for (size_t v : nodes)
		tmp.push_back(make_pair(v, r[v]));
	struct { bool operator()(pair<size_t, double> p, pair<size_t, double> q) {return p.second > q.second;} } comp;
	sort(tmp.begin(), tmp.end(), comp);
	for (size_t i = 0; i < nodes.size(); i++)
		nodes[i] = tmp[i].first;
	double sum = 0, ich = 0;
	for (size_t i = 1; i < n; i++)
	{
		sum += r[nodes[i - 1]];
		if (i == h)
			ich = 1;
		else if (i > h)
			ich = (ich * i) / (i - h);
		if (min(ich * factor, sum) - i * den >= max(1.0 / n, ceil(den * i) - (den * i)))
			return false;
	}
	return true;
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

unordered_map< size_t, vector<size_t> > residual(Graph& g, Vertex& s, Vertex& t, unordered_map<size_t, size_t>& index)
{
	boost::push_relabel_max_flow(g, s, t);
	unordered_map< size_t, vector<size_t> > adj = unordered_map< size_t, vector<size_t> >();
	for (auto it = edges(g); it.first != it.second; it.first++)
	{
		Edge e = *(it.first);
		if (get(boost::edge_residual_capacity, g)[e] > 1e-6)
		{
			size_t u = index[get(boost::vertex_index, g)[boost::source(e, g)]], v = index[get(boost::vertex_index, g)[boost::target(e, g)]];
			if (adj.find(u) == adj.end())
				adj[u] = vector<size_t>();
			adj[u].push_back(v);
		}
	}
	return adj;
}

void compute_desc(unordered_map< size_t, unordered_set<size_t> >& scc_graph, size_t u, size_t n, unordered_map< size_t, unordered_set<size_t> >& scc_nodes, unordered_map< size_t, unordered_set<size_t> >& desc, unordered_map<size_t, unordered_set<size_t> >& closed_nodes)
{
	if (desc.find(u) == desc.end())
		desc[u] = unordered_set<size_t>();
	if (closed_nodes.find(u) == closed_nodes.end())
		closed_nodes[u] = unordered_set<size_t>();
	for (size_t v : scc_nodes[u])
		if (v < n)
			closed_nodes[u].insert(v);
	for (size_t v : scc_graph[u])
	{
		if (desc.find(v) == desc.end())
			compute_desc(scc_graph, v, n, scc_nodes, desc, closed_nodes);
		desc[u].insert(desc[v].begin(), desc[v].end());
		desc[u].insert(v);
		closed_nodes[u].insert(closed_nodes[v].begin(), closed_nodes[v].end());
	}
}

void compute_scc(unordered_map< size_t, vector<size_t> >& adj, size_t n, unordered_map< size_t, unordered_set<size_t> >& scc_nodes, unordered_map< size_t, unordered_set<size_t> >& desc, unordered_map<size_t, unordered_set<size_t> >& closed_nodes)
{
	unordered_map<size_t, size_t> index = unordered_map<size_t, size_t>(), low = unordered_map<size_t, size_t>(), scc = unordered_map<size_t, size_t>();
	unordered_map<size_t, size_t> pos = unordered_map<size_t, size_t>(), parent = unordered_map<size_t, size_t>();
	stack<size_t> s = stack<size_t>();
	unordered_set<size_t> present = unordered_set<size_t>();
	size_t cnt = 0, num_scc = 0;
	for (auto e : adj)
	{
		size_t u = e.first;
		if (index.find(u) == index.end())
		{
			index[u] = cnt;
			low[u] = cnt;
			cnt++;
			pos[u] = 0;
			s.push(u);
			present.insert(u);
			size_t v = u;
			while (!s.empty())
			{
				while (pos[v] < adj[v].size())
				{
					size_t w = adj[v][pos[v]];
					if (index.find(w) == index.end())
					{
						index[w] = cnt;
						low[w] = cnt;
						cnt++;
						s.push(w);
						present.insert(w);
						parent[w] = v;
						v = w;
						break;
					}
					if (present.find(w) != present.end())
						low[v] = min(low[v], low[w]);
					pos[v]++;
				}
				if (pos[v] >= adj[v].size())
				{
					if (index[v] == low[v])
					{
						size_t w;
						do
						{
							w = s.top();
							s.pop();
							present.erase(w);
							scc[w] = num_scc;
							if (scc_nodes.find(num_scc) == scc_nodes.end())
								scc_nodes[num_scc] = unordered_set<size_t>();
							scc_nodes[num_scc].insert(w);
						}
						while (w != v);
						num_scc++;
					}
					v = parent[v];
				}
			}
		}
	}
	scc_nodes.erase(scc[n]);
	scc_nodes.erase(scc[n + 1]);
	unordered_map< size_t, unordered_set<size_t> > scc_graph = unordered_map< size_t, unordered_set<size_t> >();
	for (auto e : adj)
	{
		size_t u = e.first;
		if (scc[u] != scc[n] and scc[u] != scc[n + 1])
		{
			if (scc_graph.find(scc[u]) == scc_graph.end())
				scc_graph[scc[u]] = unordered_set<size_t>();
			for (size_t v : e.second)
				if (scc[v] != scc[u] and scc[v] != scc[n] and scc[v] != scc[n + 1])
					scc_graph[scc[u]].insert(scc[v]);
		}
	}
	adj.clear();
	for (auto c : scc_graph)
		if (desc.find(c.first) == desc.end())
			compute_desc(scc_graph, c.first, n, scc_nodes, desc, closed_nodes);
}

void enumerate_densest(unordered_set<size_t> c1, unordered_set<size_t> c2, size_t comp, size_t n, unordered_map< size_t, unordered_set<size_t> >& desc, unordered_map< size_t, unordered_set<size_t> >& anc, unordered_map< size_t, unordered_set<size_t> >& scc_nodes, unordered_map< size_t, unordered_set<size_t> >& closed_nodes, vector< vector<size_t> >& res)
{
	if (comp != numeric_limits<int>::max())
	{
		c1.insert(comp);
		for (size_t c : desc[comp])
			c2.erase(c);
		for (size_t c : anc[comp])
			c2.erase(c);
	}
	if (!c1.empty())
	{
		vector<size_t> s = vector<size_t>();
		unordered_set<size_t> v = unordered_set<size_t>();
		for (size_t c : c1)
			v.insert(closed_nodes[c].begin(), closed_nodes[c].end());
		s.insert(s.end(), v.begin(), v.end());
		res.push_back(s);
	}
	while (!c2.empty())
	{
		size_t c = *(c2.begin());
		c2.erase(c);
		bool choose = false;
		for (size_t v : scc_nodes[c])
			if (v < n)
			{
				choose = true;
				break;
			}
		if (choose)
			enumerate_densest(c1, c2, c, n, desc, anc, scc_nodes, closed_nodes, res);
	}
}

string code(vector<size_t> s)
{
	string h = "";
	sort(s.begin(), s.end());
	for (size_t v : s)
		h += (to_string(v) + '-');
	return h;
}

vector<size_t> decode(string h)
{
	vector<size_t> s = vector<size_t>();
	size_t pos = 0;
	for (size_t i = 0; i < h.size(); i++)
		if (h[i] == '-')
		{
			s.push_back(stoi(h.substr(pos, i)));
			pos = i + 1;
		}
	return s;
}
