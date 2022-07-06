#include "util.hpp"

unordered_map<size_t, size_t> core_reduce(vector< set<size_t> >& cliques, double& max_den)
{
	unordered_map<size_t, size_t> vert = unordered_map<size_t, size_t>(), deg = unordered_map<size_t, size_t>(), pos = unordered_map<size_t, size_t>();
	unordered_map< size_t, unordered_set<size_t> > adj = unordered_map< size_t, unordered_set<size_t> >();
	size_t m = cliques.size(), md = 0;
	for (size_t i = 0; i < cliques.size(); i++)
	{
		for (size_t v : cliques[i])
		{
			if (deg.find(v) == deg.end())
				deg[v] = 0;
			deg[v]++;
			if (adj.find(v) == adj.end())
				adj[v] = unordered_set<size_t>();
			adj[v].insert(i);
		}
	}
	for (auto e : deg)
		if (md < e.second)
			md = e.second;
	vector<size_t> bin = vector<size_t>(md + 1, 0);
	for (auto e : deg)
		bin[e.second]++;
	size_t start = 1, n = adj.size();
	for (size_t d = 0; d <= md; d++)
	{
		size_t num = bin[d];
		bin[d] = start;
		start += num;
	}
	for (auto d : deg)
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
		for (size_t j : adj[v])
		{
			for (size_t u : cliques[j])
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
				if (u != v)
					adj[u].erase(j);
			}
			m--;
		}
		n--;
		double den = (double) m / n;
		if (n and den > max_den)
			max_den = den;
		adj.erase(v);
	}
	return deg;
}

unordered_map<size_t, double> kclistpp(vector< set<size_t> >& cliques, size_t T, unordered_map< size_t, unordered_map<size_t, double> >& alpha)
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
		for (auto e : alpha)
			for (auto f : e.second)
				alpha[e.first][f.first] *= (1 - gamma);
		for (auto e : r)
			r[e.first] *= (1 - gamma);
		for (size_t i = 0; i < cliques.size(); i++)
		{
			size_t u = *(cliques[i].begin());
			for (size_t v : cliques[i])
				if (r[v] < r[u])
					u = v;
			alpha[u][i] += gamma;
			r[u] += gamma;
		}
	}
	return r;
}

pair<size_t, double> extract_densest(vector< set<size_t> >& cliques, vector<size_t>& nodes, unordered_map< size_t, unordered_map<size_t, double> >& alpha, unordered_map<size_t, double>& r)
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
	for (set<size_t>& c : cliques)
	{
		size_t l = 0;
		for (size_t v : c)
			if (level[v] > l)
				l = level[v];
		tentative[l]++;
	}
	pava(nodes, nag, val, tentative, level);
	for (size_t i = 0; i < cliques.size(); i++)
	{
		set<size_t> c = cliques[i];
		size_t max_level = 0, max_level_cnt = 0;
		for (size_t v : c)
			if (level[v] > max_level)
			{
				max_level = level[v];
				max_level_cnt = 1;
			}
			else if (level[v] == max_level)
				max_level_cnt++;
		double sum = 0;
		for (size_t v : c)
			if (level[v] < max_level)
			{
				sum += alpha[v][i];
				r[v] -= alpha[v][i];
				alpha[v][i] = 0;
			}
		for (size_t v : c)
			if (level[v] == max_level)
			{
				r[v] += (sum / max_level_cnt);
				alpha[v][i] += (sum / max_level_cnt);
			}
	}
	return make_pair(nag[0], val[0]);
}

bool densest_maxflow(vector< set<size_t> >& cliques, vector<size_t>& nodes, size_t n)
{
    Graph g = Graph();
    vector<Vertex> R = vector<Vertex>();
	Vertex s = boost::add_vertex(g), t = boost::add_vertex(g);
	unordered_map<size_t, size_t> ids = unordered_map<size_t, size_t>();
	for (size_t i = 0; i < nodes.size(); i++)
		ids[nodes[i]] = i;
	for (size_t i = 0; i < n; i++)
		R.push_back(boost::add_vertex(g));
	size_t m = 0;
	for (set<size_t> c : cliques)
	{
		bool flag = true;
		for (size_t v : c)
			if (ids[v] >= n)
			{
				flag = false;
				break;
			}
		if (flag)
		{
			m++;
			Vertex v = boost::add_vertex(g);
			for (size_t u : c)
				insert_edge(g, v, R[ids[u]], n);
			insert_edge(g, s, v, n);
		}
	}
	for (size_t i = 0; i < n; i++)
		insert_edge(g, R[i], t, m);
	return boost::push_relabel_max_flow(g, s, t) >= m * n * (1 - 1e-6);
}

pair<size_t, double> maximum_density(vector< set<size_t> >& cliques, size_t h, vector<size_t>& nodes)
{
	unordered_map< size_t, unordered_map<size_t, double> > alpha = unordered_map< size_t, unordered_map<size_t, double> >();
	for (size_t i = 0; i < cliques.size(); i++)
		for (size_t v : cliques[i])
		{
			if (alpha.find(v) == alpha.end())
				alpha[v] = unordered_map<size_t, double>();
			alpha[v][i] = 1.0 / h;
		}
	for (size_t t = 1; ; t *= 2)
	{
		unordered_map<size_t, double> r = kclistpp(cliques, t, alpha);
		pair<size_t, double> res = extract_densest(cliques, nodes, alpha, r);
		double prefix_min_rho = r[nodes[0]], suffix_max_rho = -1;
		for (size_t i = 1; i < res.first; ++i)
			if (prefix_min_rho > r[nodes[i]])
				prefix_min_rho = r[nodes[i]];
		for (size_t i = nodes.size() - 1; i >= res.first; --i)
			if (suffix_max_rho < r[nodes[i]])
				suffix_max_rho = r[nodes[i]];
		if (prefix_min_rho * (1 - 1e-6) > suffix_max_rho)
		{
			if (densest_goldberg(nodes, res.first, h, r, res.second, 1))
				return res;
			if (densest_maxflow(cliques, nodes, res.first))
				return res;
		}
	}
	return make_pair(0, 0);
}

unordered_map< size_t, vector<size_t> > construct(size_t source, size_t h, vector< set<size_t> >& cliques, double rho)
{
    Graph g = Graph();
	unordered_map<size_t, size_t> deg = unordered_map<size_t, size_t>(), index = unordered_map<size_t, size_t>();
	unordered_map<size_t, Vertex> rv = unordered_map<size_t, Vertex>(), rc = unordered_map<size_t, Vertex>();
	Vertex s = boost::add_vertex(g), t = boost::add_vertex(g);
	index[get(boost::vertex_index, g)[s]] = source;
	index[get(boost::vertex_index, g)[t]] = source + 1;
	for (size_t i = 0; i < cliques.size(); i++)
	{
		if (rc.find(i) == rc.end())
		{
			rc[i] = boost::add_vertex(g);
			index[get(boost::vertex_index, g)[rc[i]]] = source + 2 + i;
		}
		for (size_t v : cliques[i])
		{
			if (deg.find(v) == deg.end())
				deg[v] = 0;
			deg[v]++;
			if (rv.find(v) == rv.end())
			{
				rv[v] = boost::add_vertex(g);
				index[get(boost::vertex_index, g)[rv[v]]] = v;
			}
			insert_edge(g, rv[v], rc[i], 1, h - 1);
		}
	}
	for (auto& d : deg)
	{
		insert_edge(g, s, rv[d.first], d.second);
		insert_edge(g, rv[d.first], t, h * rho);
	}
	return residual(g, s, t, index);
}
