#include "util.hpp"

void core_reduce(unordered_map< size_t, unordered_set<size_t> >& adj, double& max_den)
{
	unordered_map<size_t, size_t> vert = unordered_map<size_t, size_t>(), deg = unordered_map<size_t, size_t>(), pos = unordered_map<size_t, size_t>();
	size_t m = 0, md = 0, n = adj.size();
	for (auto& e : adj)
	{
		m += e.second.size();
		deg[e.first] = e.second.size();
		if (md < e.second.size())
			md = e.second.size();
	}
	m /= 2;
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
	max_den = (double) m / n;
	for (size_t i = 1; i <= adj.size(); i++)
	{
		size_t v = vert[i];
		for (size_t u : adj[v])
		{
			if (pos[u] >= i)
			{
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
				m--;
			}
		}
		n--;
		double den = (double) m / n;
		if (n and den > max_den)
			max_den = den;
	}
	for (auto& d : deg)
		if (d.second >= max_den)
		{
			unordered_set<size_t> s = adj[d.first];
			for (size_t v : s)
				if (deg[v] < max_den)
					adj[d.first].erase(v);
		}
		else
			adj.erase(d.first);
}

unordered_map<size_t, double> kclistpp(unordered_map< size_t, unordered_set<size_t> >& adj, size_t T, unordered_map< size_t, unordered_map<size_t, double> >& alpha)
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
			for (size_t v : e.second)
				if (u < v)
				{
					if (r[v] < r[u])
					{
						alpha[v][u] += gamma;
						r[v] += gamma;
					}
					else
					{
						alpha[u][v] += gamma;
						r[u] += gamma;
					}
				}
		}
	}
	return r;
}

pair<size_t, double> extract_densest(unordered_map< size_t, unordered_set<size_t> >& adj, vector<size_t>& nodes, unordered_map< size_t, unordered_map<size_t, double> >& alpha, unordered_map<size_t, double>& r)
{
	vector<size_t> nag = vector<size_t>();
	vector<double> val = vector<double>();
	vector<size_t> tentative = vector<size_t>(nodes.size(), 0);
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
		for (size_t v : e.second)
			if (u < v)
			{
				if (level[v] > level[u])
					tentative[level[v]]++;
				else
					tentative[level[u]]++;
			}
	}
	pava(nodes, nag, val, tentative, level);
	for (auto& e : adj)
	{
		size_t u = e.first;
		for (size_t v : e.second)
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
	return make_pair(nag[0], val[0]);
}

bool densest_maxflow(unordered_map< size_t, unordered_set<size_t> >& adj, vector<size_t>& nodes, size_t n, double den)
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
			size_t deg = 0;
			for (size_t v : e.second)
				if (ids[v] < n)
				{
					deg++;
					if (u < v)
						insert_edge(g, R[ids[u]], R[ids[v]], 1, 1);
				}
			insert_edge(g, s, R[ids[u]], deg);
			insert_edge(g, R[ids[u]], t, 2 * den);
		}
	}
	return boost::push_relabel_max_flow(g, s, t) >= 2 * n * den * (1 - 1e-6);
}

pair<size_t, double> maximum_density(unordered_map< size_t, unordered_set<size_t> >& adj, vector<size_t>& nodes)
{
	unordered_map< size_t, unordered_map<size_t, double> > alpha = unordered_map< size_t, unordered_map<size_t, double> >();
	for (auto& e : adj)
	{
		size_t u = e.first;
		if (alpha.find(u) == alpha.end())
			alpha[u] = unordered_map<size_t, double>();
		for (size_t v : e.second)
			alpha[u][v] = 0.5;
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
		if (prefix_min_rho * (1 - 1e-6) > suffix_max_rho)
		{
			if (densest_goldberg(nodes, res.first, 2, r, res.second, 1))
				return res;
			if (densest_maxflow(adj, nodes, res.first, res.second))
				return res;
		}
	}
	return make_pair(0, 0);
}

unordered_map< size_t, vector<size_t> > construct(size_t source, unordered_map< size_t, unordered_set<size_t> >& adj, double rho)
{
    Graph g = Graph();
	unordered_map<size_t, size_t> index = unordered_map<size_t, size_t>();
	unordered_map<size_t, Vertex> R = unordered_map<size_t, Vertex>();
	Vertex s = boost::add_vertex(g), t = boost::add_vertex(g);
	index[get(boost::vertex_index, g)[s]] = source;
	index[get(boost::vertex_index, g)[t]] = source + 1;
	for (auto& e : adj)
	{
		size_t u = e.first;
		if (R.find(u) == R.end())
		{
			R[u] = boost::add_vertex(g);
			index[get(boost::vertex_index, g)[R[u]]] = u;
		}
		for (size_t v : e.second)
		{
			if (R.find(v) == R.end())
			{
				R[v] = boost::add_vertex(g);
				index[get(boost::vertex_index, g)[R[v]]] = v;
			}
			insert_edge(g, R[u], R[v], 1, 1);
		}
		insert_edge(g, s, R[u], e.second.size());
		insert_edge(g, R[u], t, 2 * rho);
	}
	adj.clear();
	return residual(g, s, t, index);
}
