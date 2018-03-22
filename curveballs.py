from random import random,sample,shuffle,randint,randrange
from numpy import array, mean, std, arange,where,zeros
from igraph import Graph
from copy import deepcopy


def edgelist_2_adjlist(edge_list,mode='d'):
	if mode != 'u':
		is_dir=True
	else:
		is_dir=False
	g=Graph.TupleList(edge_list,directed=is_dir)
	adj_list = g.get_adjlist()
	return array([set(i) for i in adj_list])


def net_from_edgelist(edge_list,mode='d'):
	if mode != 'u':
		is_dir = True
	else:
		is_dir = False
	net = Graph.TupleList(edge_list,directed=is_dir)
	return net


def net_from_adjlist(adjlist,mode='d'):
	edgelist = []
	if mode != 'u':
		is_dir = True
		for i in range(len(adjlist)):
			row = [[i,j] for j in adjlist[i]]
			edgelist += row
	else:
		is_dir = False
		for i in range(len(adjlist)):
			for j in adjlist[i]:
				if [j, i] not in edgelist:
					edgelist.append([i,j])
	return Graph.TupleList(edgelist, directed=is_dir)


def curveballs(adjlist, reps_n = 1000, mode='b',swaplike='no'):		#alternative modes are 'u': undirected; 'd': directed
	for rep in range(reps_n):
		id_a,id_b=sample(range(len(adjlist)),2)
		a, b = adjlist[id_a],adjlist[id_b]
		if mode=='b':
			ua = a - b		# unique elements of a
		else:
			ua = a - (b | {id_b})  # unique elements of a
		if len(ua) != 0:		# first check if trade is possible
			if mode == 'b':
				ub = b - a		#unique elements of b
			else:
				ub = b - (a | {id_a})  # unique elements of b
			if len(ub) != 0:		# second check if trade is possible
				if swaplike != 'no':
					trade_size = 1
				else:
					trade_size = min(len(ua),len(ub))		#allow trade size = 0?
				ra = set(sample(ua,trade_size))
				rb = set(sample(ub, trade_size))
				adjlist[id_a] = (a | rb) - ra
				adjlist[id_b] = (b | ra) - rb
				if mode == 'u':
					ra, rb = map(array,map(list,[ra,rb]))
					adjlist[rb] |= {id_a}
					adjlist[rb]	-= {id_b}
					adjlist[ra] |= {id_b}
					adjlist[ra] -= {id_a}
	return adjlist


def perturbation_score_edge_list(e_r,e_0,mode='d'):
	if mode == 'u':
		e_r = set(['_'.join(map(str, i)) for i in e_r]+['_'.join(map(str, i[::-1])) for i in e_r])
		e_0 = set(['_'.join(map(str, i)) for i in e_0]+['_'.join(map(str, i[::-1])) for i in e_0])
	else:
		e_r = set(['_'.join(map(str, i)) for i in e_r])
		e_0 = set(['_'.join(map(str, i)) for i in e_0])
	return 1-len(e_r & e_0)/float(len(e_0))


def perturbation_score_adj_list(al_r, al_0):
	return sum([len(al_r[i]-al_0[i]) for i in range(len(al_0))])/float(sum([len(j) for j in al_0]))


#check that degree is preserved
for test in range(100):
	r_mode=sample(['d','u','b'],1)[0]
	if r_mode == 'd':
		ee_net = []
		while len(ee_net) < 1000:
			edge = sample(range(100), 2)
			if edge not in ee_net:
				ee_net.append(edge)
	elif r_mode == 'b':
		ee_net = []
		while len(ee_net)<1000:
			edge = [randrange(50),randrange(50,100)]
			if edge not in ee_net:
				ee_net.append(edge)
	else:
		ee_net = []
		while len(ee_net)<1000:
			edge = sample(range(100),2)
			if edge[::-1] not in ee_net and edge not in ee_net:
				ee_net.append(edge)
	net_0 = net_from_edgelist(ee_net,mode=r_mode)
	al = edgelist_2_adjlist(ee_net,mode=r_mode)
	al_0 = deepcopy(al)
	al_r = curveballs(al,reps_n=100,mode=r_mode)
	#print perturbation_score_adj_list(al_r,al_0)
	net_r = net_from_adjlist(al_r,mode=r_mode)
	if net_0 != net_r and sorted(net_0.degree(mode = 'IN')) == sorted(net_r.degree(mode = 'IN')) and sorted(net_0.degree(mode = 'OUT')) == sorted(net_r.degree(mode = 'OUT')):
		print test,'OK',r_mode
	else:
		break


#perturbation experiments
def make_er_net(max_nodes,edge_n,mode='d'):
	r_mode = 'd'
	if r_mode == 'd':
		er_net = []
		while len(er_net) < 10000:
			edge = sample(range(1000), 2)
			if edge not in er_net:
				er_net.append(edge)
	elif r_mode == 'b':
		er_net = []
		while len(er_net) < 10000:
			edge = [randrange(500), randrange(500, 1000)]
			if edge not in er_net:
				er_net.append(edge)
	else:
		er_net = []
		while len(er_net) < 10000:
			edge = sample(range(1000), 2)
			if edge[::-1] not in er_net and edge not in er_net:
				er_net.append(edge)
	return er_net


#er networks
for r_mode in ['d', 'u']:
	out = open(r_mode + '_ee.csv', 'w')
	for rep in range(100):
		nodes = randrange(100,1000)
		edges = randrange(5,50)*nodes
		er_net = make_er_net(nodes, edges, mode=r_mode)
		net_0 = deepcopy(er_net)
		al = edgelist_2_adjlist(er_net, mode=r_mode)
		al_0 = deepcopy(al)
		al_swap = deepcopy(al)
		out.write(','.join(map(str, [0, 0.0, 0.0])) + '\n')
		for step in range(1,250):
			al = curveballs(al, reps_n=100, mode=r_mode)
			al_swap = curveballs(al_swap, reps_n=100, mode=r_mode, swaplike='yes')
			sc_trade = perturbation_score_adj_list(al, al_0)
			sc_swap = perturbation_score_adj_list(al_swap, al_0)
			out.write(','.join(map(str,[step*100, sc_trade, sc_swap])) + '\n')
		print r_mode, rep
	out.close()



#barabasi

for r_mode in ['d', 'u']:
	out = open(r_mode + '_ba.csv', 'w')
	for rep in range(100):
		nodes = randrange(100,1000)
		g = Graph.Barabasi(nodes)
		ba_net = [list(i.tuple) for i in g.es]
		net_0 = deepcopy(ba_net)
		al = edgelist_2_adjlist(ba_net, mode=r_mode)
		al_0 = deepcopy(al)
		al_swap = deepcopy(al)
		out.write(','.join(map(str, [0, 0.0, 0.0])) + '\n')
		for step in range(1,250):
			al = curveballs(al, reps_n=100, mode=r_mode)
			al_swap = curveballs(al_swap, reps_n=100, mode=r_mode, swaplike='yes')
			sc_trade = perturbation_score_adj_list(al, al_0)
			sc_swap = perturbation_score_adj_list(al_swap, al_0)
			out.write(','.join(map(str,[step*100, sc_trade, sc_swap])) + '\n')
		print r_mode, rep
	out.close()




#real networks
#undirected co-occurrence network

cooc_net_file = open('cooc.csv','r')
cooc_edge_list = [i.split(',') for i in cooc_net_file]
cooc_net_file.close()
cooc_net_g = Graph.TupleList(cooc_edge_list, directed = False)
cooc_net = [i.tuple for i in cooc_net_g.es]
r_mode = 'u'
out = open('cooc_net_res.csv', 'w')
for rep in range(100):
	net_0 = deepcopy(cooc_net)
	al = edgelist_2_adjlist(cooc_net, mode=r_mode)
	al_0 = deepcopy(al)
	al_swap = deepcopy(al)
	out.write(','.join(map(str, [0, 0.0, 0.0])) + '\n')
	for step in range(1,100):
		al = curveballs(al, reps_n=100, mode=r_mode)
		al_swap = curveballs(al_swap, reps_n=100, mode=r_mode, swaplike='yes')
		sc_trade = perturbation_score_adj_list(al, al_0)
		sc_swap = perturbation_score_adj_list(al_swap, al_0)
		out.write(','.join(map(str,[step*100, sc_trade, sc_swap])) + '\n')
	print r_mode, rep


out.close()



#food web
food_web_file = open('food_web.csv','r')
food_web_edge_list = [i.split(',') for i in food_web_file]
food_web_file.close()
food_web_g = Graph.TupleList(food_web_edge_list, directed = False)
food_web = [i.tuple for i in food_web_g.es]
r_mode = 'd'
out = open('food_web_res.csv', 'w')
for rep in range(100):
	net_0 = deepcopy(food_web)
	al = edgelist_2_adjlist(food_web, mode=r_mode)
	al_0 = deepcopy(al)
	al_swap = deepcopy(al)
	out.write(','.join(map(str, [0, 0.0, 0.0])) + '\n')
	for step in range(1,250):
		al = curveballs(al, reps_n=100, mode=r_mode)
		al_swap = curveballs(al_swap, reps_n=100, mode=r_mode, swaplike='yes')
		sc_trade = perturbation_score_adj_list(al, al_0)
		sc_swap = perturbation_score_adj_list(al_swap, al_0)
		out.write(','.join(map(str,[step*100, sc_trade, sc_swap])) + '\n')
	print r_mode, rep


out.close()







