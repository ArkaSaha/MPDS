#include "util.hpp"

unordered_map<string, size_t> get_patterns(unordered_map< size_t, unordered_set<size_t> >& graph, string pattern)
{
	unordered_map<string, size_t> instances = unordered_map<string, size_t>();
	if (pattern == "diamond")
	{
		for (auto& e : graph)
		{
			size_t v = e.first;
			unordered_map< size_t, vector<size_t> > paths = unordered_map< size_t, vector<size_t> >();
			for (size_t u : e.second)
				for (size_t w : graph[u])
					if (w != v)
					{
						if (paths.find(w) == paths.end())
							paths[w] = vector<size_t>();
						paths[w].push_back(u);
					}
			for (auto& p : paths)
			{
				size_t u = p.first;
				for (size_t i = 0; i < p.second.size(); i++)
					for (size_t j = i + 1; j < p.second.size(); j++)
					{
						string str = code(vector<size_t>({v, u, p.second[i], p.second[j]}));
						if (instances.find(str) == instances.end())
							instances[str] = 0;
						instances[str]++;
					}
			}
		}
		for (auto& i : instances)
			instances[i.first] /= 4;
	}
	else if (pattern == "2-star")
		for (auto& e : graph)
		{
			size_t v = e.first;
			vector<size_t> nbr = vector<size_t>(e.second.begin(), e.second.end());
			for (size_t i = 0; i < e.second.size(); i++)
				for (size_t j = i + 1; j < e.second.size(); j++)
				{
					string str = code(vector<size_t>({v, nbr[i], nbr[j]}));
					if (instances.find(str) == instances.end())
						instances[str] = 0;
					instances[str]++;
				}
		}
	else if (pattern == "3-star")
		for (auto& e : graph)
		{
			size_t v = e.first;
			vector<size_t> nbr = vector<size_t>(e.second.begin(), e.second.end());
			for (size_t i = 0; i < e.second.size(); i++)
				for (size_t j = i + 1; j < e.second.size(); j++)
					for (size_t k = j + 1; k < e.second.size(); k++)
					{
						string str = code(vector<size_t>({v, nbr[i], nbr[j], nbr[k]}));
						if (instances.find(str) == instances.end())
							instances[str] = 0;
						instances[str]++;
					}
		}
	else if (pattern == "c3-star")
	{
		vector< set<size_t> > cliques = get_cliques(graph, 3);
		for (set<size_t> c : cliques)
			for (size_t v : c)
				for (size_t u : graph[v])
					if (c.find(u) == c.end())
					{
						vector<size_t> w = vector<size_t>(c.begin(), c.end());
						w.push_back(u);
						string str = code(w);
						if (instances.find(str) == instances.end())
							instances[str] = 0;
						instances[str]++;
					}
	}
	else
		assert(false);
	return instances;
}

void core_reduce_diamond(unordered_map< size_t, unordered_set<size_t> >& graph, vector< vector<size_t> >& choose)
{
	unordered_map<size_t, size_t> pat_deg = unordered_map<size_t, size_t>(), core = unordered_map<size_t, size_t>();
	size_t m = 0, n = graph.size();
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
	while (!heap.empty())
	{
		node nd = heap.top();
		heap.pop();
		size_t v = nd.vertex, d = nd.degree;
		handles.erase(v);
		core[v] = d;
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
		if (n and den > max_den)
			max_den = den;
		pat_deg.erase(v);
	}
	for (auto& d : core)
	{
		size_t v = d.first;
		if (d.second >= max_den)
		{
			unordered_set<size_t> nbr = graph[v];
			for (size_t u : nbr)
				if (core[u] < max_den)
					graph[v].erase(u);
		}
		else
			graph.erase(v);
	}
}

void core_reduce_star(unordered_map< size_t, unordered_set<size_t> >& graph, vector< vector<size_t> >& choose, size_t x)
{
	unordered_map<size_t, size_t> deg = unordered_map<size_t, size_t>(), pat_deg = unordered_map<size_t, size_t>(), core = unordered_map<size_t, size_t>();
	size_t m = 0, n = graph.size();
	for (auto& e : graph)
	{
		size_t v = e.first;
		deg[v] = e.second.size();
		pat_deg[v] = choose[deg[v]][x];
		for (size_t u : e.second)
			if (!graph[u].empty())
				pat_deg[v] += choose[graph[u].size() - 1][x - 1];
		m += choose[deg[v]][x];
	}
	double max_den = (double) m / n;
	fib_heap heap = fib_heap();
	unordered_map<size_t, handle_t> handles = unordered_map<size_t, handle_t>();
	for (auto& d : pat_deg)
		handles[d.first] = heap.push(node(d.first, d.second));
	while (!heap.empty())
	{
		node nd = heap.top();
		heap.pop();
		size_t v = nd.vertex, d = nd.degree;
		handles.erase(v);
		core[v] = d;
		size_t count = choose[deg[v]][x];
		for (size_t u : graph[v])
			if (handles.find(u) != handles.end() and deg[u])
				count += choose[deg[u] - 1][x - 1];
		for (size_t u : graph[v])
			if (handles.find(u) != handles.end() and deg[u] and deg[v])
			{
				heap.update(handles[u], node(u, max((*handles[u]).degree - (choose[deg[u] - 1][x - 1] + choose[deg[v] - 1][x - 1]), d)));
				for (size_t w : graph[u])
					if (handles.find(w) != handles.end() and w != v and deg[u] >= 2)
						heap.update(handles[w], node(w, max((*handles[w]).degree - choose[deg[u] - 2][x - 2], d)));
			}
		m -= count;
		n--;
		double den = (double) m / n;
		if (n and den > max_den)
			max_den = den;
		for (size_t u : graph[v])
			deg[u]--;
		deg.erase(v);
		pat_deg.erase(v);
	}
	for (auto& d : core)
	{
		size_t v = d.first;
		if (d.second >= max_den)
		{
			unordered_set<size_t> nbr = graph[v];
			for (size_t u : nbr)
				if (core[u] < max_den)
					graph[v].erase(u);
		}
		else
			graph.erase(v);
	}
}

void core_reduce_cstar(unordered_map< size_t, unordered_set<size_t> >& graph)
{
	unordered_map<size_t, size_t> pat_deg = unordered_map<size_t, size_t>(), core = unordered_map<size_t, size_t>();
	unordered_map< size_t, unordered_set<size_t> > adj = unordered_map< size_t, unordered_set<size_t> >();
	size_t m = 0;
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
	size_t n = pat_deg.size();
	double max_den = (double) m / n;
	fib_heap heap = fib_heap();
	unordered_map<size_t, handle_t> handles = unordered_map<size_t, handle_t>();
	for (auto& d : pat_deg)
		handles[d.first] = heap.push(node(d.first, d.second));
	while (!heap.empty())
	{
		node nd = heap.top();
		heap.pop();
		size_t v = nd.vertex, d = nd.degree;
		handles.erase(v);
		core[v] = d;
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
		if (n and den > max_den)
			max_den = den;
		pat_deg.erase(v);
		adj.erase(v);
	}
	vector<size_t> rem = vector<size_t>();
	for (auto& e : graph)
		if (core.find(e.first) == core.end())
			rem.push_back(e.first);
	while (!rem.empty())
	{
		size_t v = rem.back();
		rem.pop_back();
		graph.erase(v);
	}
	for (auto& d : core)
	{
		size_t v = d.first;
		if (d.second >= max_den)
		{
			unordered_set<size_t> nbr = graph[v];
			for (size_t u : nbr)
				if (core[u] < max_den)
					graph[v].erase(u);
		}
		else
			graph.erase(v);
	}
}

unordered_map<size_t, double> kclistpp(unordered_map<string, size_t>& instances, size_t T, unordered_map< size_t, unordered_map<string, double> >& alpha)
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
		for (auto& i : instances)
		{
			vector<size_t> p = decode(i.first);
			size_t u = p[0];
			for (size_t v : p)
				if (r[v] < r[u])
					u = v;
			alpha[u][i.first] += (gamma * i.second);
			r[u] += (gamma * i.second);
		}
	}
	return r;
}

pair<size_t, double> extract_densest(unordered_map<string, size_t>& instances, vector<size_t>& nodes, unordered_map< size_t, unordered_map<string, double> >& alpha, unordered_map<size_t, double>& r)
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
	for (auto& i : instances)
	{
		vector<size_t> p = decode(i.first);
		size_t l = 0;
		for (size_t v : p)
			if (level[v] > l)
				l = level[v];
		tentative[l] += i.second;
	}
	pava(nodes, nag, val, tentative, level);
	for (auto& i : instances)
	{
		vector<size_t> p = decode(i.first);
		size_t max_level = 0, max_level_cnt = 0;
		for (size_t v : p)
			if (level[v] > max_level)
			{
				max_level = level[v];
				max_level_cnt = 1;
			}
			else if (level[v] == max_level)
				max_level_cnt++;
		double sum = 0;
		for (size_t v : p)
			if (level[v] < max_level)
			{
				sum += alpha[v][i.first];
				r[v] -= alpha[v][i.first];
				alpha[v][i.first] = 0;
			}
		for (size_t v : p)
			if (level[v] == max_level)
			{
				r[v] += (sum / max_level_cnt);
				alpha[v][i.first] += (sum / max_level_cnt);
			}
	}
	return make_pair(nag[0], val[0]);
}

bool densest_maxflow(unordered_map<string, size_t>& instances, vector<size_t>& nodes, size_t n)
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
	for (auto& i : instances)
	{
		vector<size_t> p = decode(i.first);
		bool flag = true;
		for (size_t v : p)
			if (ids[v] >= n)
			{
				flag = false;
				break;
			}
		if (flag)
		{
			m += i.second;
			Vertex v = boost::add_vertex(g);
			for (size_t u : p)
				insert_edge(g, v, R[ids[u]], n * i.second);
			insert_edge(g, s, v, n * i.second);
		}
	}
	for (size_t i = 0; i < n; i++)
		insert_edge(g, R[i], t, m);
	return boost::push_relabel_max_flow(g, s, t) >= m * n * (1 - 1e-6);
}

pair<size_t, double> maximum_density(unordered_map<string, size_t>& instances, size_t h, vector<size_t>& nodes, size_t factor)
{
	unordered_map< size_t, unordered_map<string, double> > alpha = unordered_map< size_t, unordered_map<string, double> >();
	for (auto& i : instances)
	{
		vector<size_t> p = decode(i.first);
		for (size_t v : p)
		{
			if (alpha.find(v) == alpha.end())
				alpha[v] = unordered_map<string, double>();
			alpha[v][i.first] = (double) i.second / h;
		}
	}
	for (size_t t = 1; ; t *= 2)
	{
		unordered_map<size_t, double> r = kclistpp(instances, t, alpha);
		pair<size_t, double> res = extract_densest(instances, nodes, alpha, r);
		double prefix_min_rho = r[nodes[0]], suffix_max_rho = -1;
		for (size_t i = 1; i < res.first; ++i)
			if (prefix_min_rho > r[nodes[i]])
				prefix_min_rho = r[nodes[i]];
		for (size_t i = nodes.size() - 1; i >= res.first; --i)
			if (suffix_max_rho < r[nodes[i]])
				suffix_max_rho = r[nodes[i]];
		if (prefix_min_rho * (1 - 1e-6) > suffix_max_rho)
		{
			if (densest_goldberg(nodes, res.first, h, r, res.second, factor))
				return res;
			if (densest_maxflow(instances, nodes, res.first))
				return res;
		}
	}
	return make_pair(0, 0);
}

unordered_map< size_t, vector<size_t> > construct(size_t source, size_t h, unordered_map<string, size_t>& instances, double rho)
{
    Graph g = Graph();
	unordered_map<size_t, size_t> deg = unordered_map<size_t, size_t>(), index = unordered_map<size_t, size_t>();
	unordered_map<size_t, Vertex> rv = unordered_map<size_t, Vertex>(), rp = unordered_map<size_t, Vertex>();
	Vertex s = boost::add_vertex(g), t = boost::add_vertex(g);
	index[get(boost::vertex_index, g)[s]] = source;
	index[get(boost::vertex_index, g)[t]] = source + 1;
	size_t num = 0;
	for (auto& i : instances)
	{
		if (rp.find(num) == rp.end())
		{
			rp[num] = boost::add_vertex(g);
			index[get(boost::vertex_index, g)[rp[num]]] = source + 2 + num;
		}
		vector<size_t> p = decode(i.first);
		for (size_t v : p)
		{
			if (deg.find(v) == deg.end())
				deg[v] = 0;
			deg[v] += i.second;
			if (rv.find(v) == rv.end())
			{
				rv[v] = boost::add_vertex(g);
				index[get(boost::vertex_index, g)[rv[v]]] = v;
			}
			insert_edge(g, rv[v], rp[num], i.second, (h - 1) * i.second);
		}
		num++;
	}
	for (auto& d : deg)
	{
		insert_edge(g, s, rv[d.first], d.second);
		insert_edge(g, rv[d.first], t, h * rho);
	}
	return residual(g, s, t, index);
}
