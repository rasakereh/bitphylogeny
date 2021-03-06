import sys
import scipy.stats
import copy

from time         import *
from numpy        import *
from numpy.random import *
from igraph       import *

from bitphylogeny.util import *

class TSSB(object):

    min_dp_alpha    = 0.001
    max_dp_alpha    = 10.0
    min_dp_gamma    = 0.001
    max_dp_gamma    = 10.0
    min_alpha_decay = 0.001
    max_alpha_decay = 0.80
    
    def __init__(self, dp_alpha=1.0, dp_gamma=1.0, root_node=None, data=None,
                 min_depth=0, max_depth=15, alpha_decay=1.0):
        if root_node is None:
            raise Exception("Root node must be specified.")
        
        self.min_depth   = min_depth
        self.max_depth   = max_depth
        self.dp_alpha    = dp_alpha
        self.dp_gamma    = dp_gamma
        self.alpha_decay = alpha_decay
        self.data        = data
        self.num_data    = 0 if data is None else data.shape[0]
        self.root        = { 'node'     : root_node,
                             'main'     : boundbeta(1.0, dp_alpha) if self.min_depth == 0 else 0.0,
                             'sticks'   : empty((0,1)),
                             'children' : [] }
        root_node.tssb = self

        if False:
            data_u           = rand(self.num_data)
            self.assignments = []
            for n in range(self.num_data):
                (c, path) = self.find_node(data_u[n])
                c.add_datum(n)
                self.assignments.append(c)
        else:
            self.assignments = []
            for n in range(self.num_data):
                self.root['node'].add_datum(n)
                self.assignments.append(self.root['node'])

    def add_data(self, data):
        (weights, nodes) = self.get_mixture()
        num_new_data = data.shape[0]
        for n in range(num_new_data):
            logprobs = []
            for k, node in enumerate(nodes):
                logprobs.append( log(weights[k]) + node.logprob(data[n]) )
            logprobs = array(logprobs)
            probs    = exp(logprobs - logsumexp(logprobs))
            best_k   = sum( rand() > cumsum(probs) )
            nodes[best_k].add_datum(n + self.num_data)
            self.assignments.append(nodes[best_k])
        self.data = vstack([self.data, data])
        self.num_data += num_new_data

    def clear_data(self):
        dims = self.data.shape[1]
        for n in range(self.num_data):
            self.assignments[n].remove_datum(n)
        self.assignments = []
        self.data        = empty((0,dims))
        self.num_data    = 0

    def resample_node_params(self, iters=1):
        for iter in range(iters):
            def descend(root):
                for index, child in enumerate(root['children']):
                    descend(child)
                root['node'].resample_params()
            descend(self.root)

    def resample_assignments(self):

        def path_lt(path1, path2):
            if len(path1) == 0 and len(path2) == 0:
                return 0
            elif len(path1) == 0:
                return 1
            elif len(path2) == 0:
                return -1
            s1 = "".join(map(lambda i: "%03d" % (i), path1))
            s2 = "".join(map(lambda i: "%03d" % (i), path2))

            return cmp(s2, s1)

        epsilon = finfo(float64).eps
        lengths = []
        reassign = 0
        better = 0        
        for n in range(self.num_data):

            # Get an initial uniform variate.
            ancestors = self.assignments[n].get_ancestors()
            current   = self.root
            indices   = []
            for anc in ancestors[1:]:
                index     = map( lambda c: c['node'], current['children']).index(anc)
                current   = current['children'][index]
                indices.append(index)
            
            max_u = 1.0
            min_u = 0.0
            llh_s = log(rand()) + self.assignments[n].logprob(self.data[n:n+1])
            #llh_s = self.assignments[n].logprob(self.data[n:n+1]) - 0.0000001
            while True:
                new_u                = (max_u-min_u)*rand() + min_u
                (new_node, new_path) = self.find_node(new_u)
                new_llh              = new_node.logprob(self.data[n:n+1])
                if new_llh > llh_s:
                    if new_node != self.assignments[n]:
                        if (new_llh > self.assignments[n].logprob(self.data[n:n+1]) ):
                            better += 1
                        self.assignments[n].remove_datum(n)
                        new_node.add_datum(n)
                        self.assignments[n] = new_node
                        reassign += 1
                    break
                elif abs(max_u-min_u) < epsilon:
                    print >>sys.stderr, "Slice sampler shrank down.  Keep current state."
                    break
                else:
                    path_comp = path_lt(indices, new_path)
                    if path_comp < 0:
                        min_u = new_u
                    elif path_comp > 0:
                        max_u = new_u
                    else:
                        raise Exception("Slice sampler weirdness.")
            lengths.append(len(new_path))
        lengths = array(lengths)
        #print "reassign: "+str(reassign)+" better: "+str(better)

    def cull_tree(self):
        def descend(root):
            counts = array(map(lambda child: descend(child), root['children']))
            keep   = len(trim_zeros(counts, 'b'))

            for child in root['children'][keep:]:
                child['node'].kill()
                del child['node']

            root['sticks']   = root['sticks'][:keep]
            root['children'] = root['children'][:keep]

            return sum(counts) + root['node'].num_local_data()

        descend(self.root)

    def resample_sticks(self):
        def descend(root, depth=0):

            data_down = 0
            indices   = range(len(root['children']))
            indices.reverse()
            for i in indices:
                child             = root['children'][i]
                child_data        = descend(child, depth+1)
                post_alpha        = 1.0 + child_data
                post_beta         = self.dp_gamma + data_down
                root['sticks'][i] = boundbeta(post_alpha, post_beta)
                data_down += child_data

            # Resample the main break.
            data_here    = root['node'].num_local_data()
            post_alpha   = 1.0 + data_here
            post_beta    = (self.alpha_decay**depth)*self.dp_alpha + data_down
            root['main'] = boundbeta( post_alpha, post_beta ) if self.min_depth <= depth else 0.0

            return data_here + data_down

        descend(self.root)

    def resample_stick_orders(self):
        def descend(root, depth=0):
            if not root['children']:
                return
           
            new_order   = []
            represented = set(filter(lambda i: root['children'][i]['node'].has_data(), 
                                     range(len(root['children']))))
            all_weights = diff(hstack([0.0, sticks_to_edges(root['sticks'])]))
            while True:
                if not represented:
                    break

                u = rand()
                while True:
                    sub_indices = filter(lambda i: i not in new_order, range(root['sticks'].shape[0]))
                    sub_weights = hstack([all_weights[sub_indices], 1.0 - sum(all_weights)])
                    sub_weights = sub_weights / sum(sub_weights)
                    index       = sum(u > cumsum(sub_weights))

                    if index == len(sub_indices):
                        root['sticks'] = vstack([ root['sticks'], boundbeta(1, self.dp_gamma) ])
                        root['children'].append({ 'node'     : root['node'].spawn(), 
                                                  'main'     : boundbeta(1.0, (self.alpha_decay**(depth+1))*self.dp_alpha) if self.min_depth <= (depth+1) else 0.0,
                                                  'sticks'   : empty((0,1)),
                                                  'children' : [] })
                        all_weights = diff(hstack([0.0, sticks_to_edges(root['sticks'])]))
                    else:
                        index = sub_indices[index]
                        break
                new_order.append(index)
                represented.discard(index)

            new_children = []
            for k in new_order:
                child = root['children'][k]
                new_children.append(child)
                descend(child, depth + 1)

            for k in filter(lambda k: k not in new_order, range(root['sticks'].shape[0])):
                root['children'][k]['node'].kill()
                del root['children'][k]['node']

            root['children'] = new_children
            root['sticks']   = zeros((len(root['children']),1))
        descend(self.root)
        
        # Immediately resample sticks.
        self.resample_sticks()

    def resample_hypers(self, dp_alpha=True, alpha_decay=True, dp_gamma=True):

        def dp_alpha_llh(dp_alpha, alpha_decay):
            def descend(dp_alpha, root, depth=0):
                llh = betapdfln(root['main'], 1.0, (alpha_decay**depth)*dp_alpha) if self.min_depth <= depth else 0.0
                for child in root['children']:
                    llh += descend(dp_alpha, child, depth+1)
                return llh            
            return descend(dp_alpha, self.root)

        if dp_alpha:
            upper = self.max_dp_alpha
            lower = self.min_dp_alpha
            llh_s = log(rand()) + dp_alpha_llh(self.dp_alpha, self.alpha_decay)
            while True:
                new_dp_alpha = (upper-lower)*rand() + lower
                new_llh       = dp_alpha_llh(new_dp_alpha, self.alpha_decay)
                if new_llh > llh_s:
                    break
                elif new_dp_alpha < self.dp_alpha:
                    lower = new_dp_alpha
                elif new_dp_alpha > self.dp_alpha:
                    upper = new_dp_alpha
                else:
                    raise Exception("Slice sampler shrank to zero!")
            self.dp_alpha = new_dp_alpha

        if alpha_decay:
            upper = self.max_alpha_decay
            lower = self.min_alpha_decay
            llh_s = log(rand()) + dp_alpha_llh(self.dp_alpha, self.alpha_decay)
            while True:
                new_alpha_decay = (upper-lower)*rand() + lower
                new_llh         = dp_alpha_llh(self.dp_alpha, new_alpha_decay)
                if new_llh > llh_s:
                    break
                elif new_alpha_decay < self.alpha_decay:
                    lower = new_alpha_decay
                elif new_alpha_decay > self.alpha_decay:
                    upper = new_alpha_decay
                else:
                    raise Exception("Slice sampler shrank to zero!")
            self.alpha_decay = new_alpha_decay

        def dp_gamma_llh(dp_gamma):
            def descend(dp_gamma, root):
                llh = 0
                for i, child in enumerate(root['children']):
                    llh += betapdfln(root['sticks'][i], 1.0, dp_gamma)
                    llh += descend(dp_gamma, child)
                return llh
            return descend(dp_gamma, self.root)

        if dp_gamma:
            upper = self.max_dp_gamma
            lower = self.min_dp_gamma
            llh_s = log(rand()) + dp_gamma_llh(self.dp_gamma)
            while True:
                new_dp_gamma = (upper-lower)*rand() + lower
                new_llh       = dp_gamma_llh(new_dp_gamma)
                if new_llh > llh_s:
                    break
                elif new_dp_gamma < self.dp_gamma:
                    lower = new_dp_gamma
                elif new_dp_gamma > self.dp_gamma:
                    upper = new_dp_gamma
                else:
                    raise Exception("Slice sampler shrank to zero!")
            self.dp_gamma = new_dp_gamma
        
    def draw_data(self, num_data=1, **args):
        self.data        = []
        self.assignments = []
        for n in range(num_data):
            u    = rand()
            (node, path) = self.find_node(u)
            self.data.append(node.sample(args))
            self.assignments.append(node)
            node.add_datum(n)
            self.num_data += 1
        self.data = concatenate(self.data)
        return self.data

    def resample_data(self, **args):
        for n in range(self.num_data):
            u    = rand()
            (node, path) = self.find_node(u)
            self.assignments[n].remove_datum(n)
            node.add_datum(n)
            self.assignments[n] = node
            self.data[n] = node.sample(args)[0]

    def find_node(self, u):
        def descend(root, u, depth=0):
            if depth >= self.max_depth:
                #print >>sys.stderr, "WARNING: Reached maximum depth."
                return (root['node'], [])
            elif u < root['main']:
                return (root['node'], [])
            else:
                # Rescale the uniform variate to the remaining interval.
                u = (u - root['main']) / (1.0 - root['main'])
                
                # Perhaps break sticks out appropriately.
                while not root['children'] or (1.0 - prod(1.0 - root['sticks'])) < u:
                    root['sticks'] = vstack([ root['sticks'], boundbeta(1, self.dp_gamma) ])
                    root['children'].append({ 'node'     : root['node'].spawn(),
                                              'main'     : boundbeta(1.0, (self.alpha_decay**(depth+1))*self.dp_alpha) if self.min_depth <= (depth+1) else 0.0,
                                              'sticks'   : empty((0,1)),
                                              'children' : [] })

                edges = 1.0 - cumprod(1.0 - root['sticks'])
                index = sum(u > edges)
                edges = hstack([0.0, edges])
                u     = (u - edges[index]) / (edges[index+1] - edges[index])

                (node, path) = descend(root['children'][index], u, depth+1)
                
                path.insert(0, index)
                
                return (node, path)
        return descend(self.root, u)

    def get_mixture(self):
        def descend(root, mass):
            weight  = [ mass * root['main'] ]
            node    = [ root['node'] ]
            edges   = sticks_to_edges(root['sticks'])
            weights = diff(hstack([0.0, edges]))

            for i, child in enumerate(root['children']):
                (child_weights, child_nodes) = descend(child, mass*(1.0-root['main'])*weights[i])
                weight.extend(child_weights)
                node.extend(child_nodes)
            return (weight, node)
        return descend(self.root, 1.0)

    def complete_data_log_likelihood(self):
        weights, nodes = self.get_mixture();
        llhs = []
        for i, node in enumerate(nodes):
            if node.num_local_data():
                llhs.append(node.num_local_data()*log(weights[i]) + node.data_log_likelihood())
        return sum(array(llhs))
    
    def complete_data_log_likelihood_nomix(self):
        weights, nodes = self.get_mixture();
        llhs = []
        lln = []
        for i, node in enumerate(nodes):
            if node.num_local_data():
                llhs.append(node.data_log_likelihood())
                lln.append(node.num_local_data()*log(weights[i]))
                #lln.append(weights[i])
        return (sum(array(lln)),sum(array(llhs)))
    
    def unnormalized_postertior_with_hypers(self):
        weights, nodes = self.get_mixture();
        llhs = []
        for i, node in enumerate(nodes):
            if node.num_local_data():
                llhs.append(node.num_local_data()*log(weights[i]) +
                            node.data_log_likelihood() +
                            node.parameter_log_prior() ) 
        llh = sum(array(llhs))
        
        def alpha_descend(root, depth=0):
            llh = betapdfln(root['main'], 1.0,
                            (self.alpha_decay**depth)*self.dp_alpha) if self.min_depth <= depth else 0.0
            for child in root['children']:
                llh += alpha_descend(child, depth+1)
            return llh            
        weights_log_prob = alpha_descend(self.root)
        
        def gamma_descend(root):
            llh = 0
            for i, child in enumerate(root['children']):
                llh += betapdfln(root['sticks'][i], 1.0, self.dp_gamma)
                llh += gamma_descend(child)
            return llh
        sticks_log_prob =  gamma_descend(self.root)
        
        return llh + weights_log_prob + sticks_log_prob
    
    def unnormalized_postertior(self):
        weights, nodes = self.get_mixture();
        llhs = []
        for i, node in enumerate(nodes):
            if node.num_local_data():
                llhs.append(node.num_local_data()*log(weights[i]) +
                            node.data_log_likelihood() +
                            node.parameter_log_prior() ) 
        llh = sum(array(llhs))
        
        return llh
    
    def resample_tree_topology(self):
        #x = self.complete_data_log_likelihood_nomix()
        post = log(rand()) + self.unnormalized_postertior()
        weights, nodes = self.get_mixture();
        if len(nodes)>1:
            nodeAnum = randint(0,len(nodes))
            nodeBnum = randint(0,len(nodes))
            
            while nodeAnum == nodeBnum:
                nodeBnum = randint(0, len(nodes))
            
            def swap_nodes(nodeAnum, nodeBnum):
                def findNodes(root, nodeNum, nodeA=False, nodeB=False):
                    node = root
                    if nodeNum == nodeAnum:
                        nodeA = node
                    if nodeNum == nodeBnum:
                        nodeB = node
                    for i, child in enumerate(root['children']):
                        nodeNum = nodeNum + 1
                        (nodeA, nodeB, nodeNum) = findNodes(child, nodeNum, nodeA, nodeB)
                    return (nodeA, nodeB, nodeNum)
                
                (nodeA,nodeB,nodeNum) = findNodes(self.root, nodeNum=0)
                
                paramsA = nodeA['node'].params
                dataA = set(nodeA['node'].data)
                mainA = nodeA['main']
            
                nodeA['node'].params = nodeB['node'].params
                
                for dataid in list(dataA):
                    nodeA['node'].remove_datum(dataid)
                for dataid in nodeB['node'].data:
                    nodeA['node'].add_datum(dataid)
                    self.assignments[dataid] = nodeA['node']
                nodeA['main']=nodeB['main']
                
                nodeB['node'].params = paramsA
                dataB = set(nodeB['node'].data)
                
                for dataid in dataB:
                    nodeB['node'].remove_datum(dataid)
                for dataid in dataA:
                    nodeB['node'].add_datum(dataid)
                    self.assignments[dataid] = nodeB['node']
                nodeB['main']=mainA
                
            
            swap_nodes(nodeAnum,nodeBnum)
            self.resample_sticks()
            post_new = self.unnormalized_postertior()
            #xn = self.complete_data_log_likelihood_nomix()
            if(post_new < post):
                swap_nodes(nodeAnum,nodeBnum)
                self.resample_sticks()
            else:
                print "successful swap!!!"


    def resample_tree_topology_root(self):
        #x = self.complete_data_log_likelihood_nomix()
        post = log(rand()) + self.unnormalized_postertior()
        weights, nodes = self.get_mixture();

        empty_root = False
        if len(nodes)>1:
            if len(nodes[0].data) == 0:
                print 'swapping root'
                empty_root = True
                nodeAnum = 0
            else:
                nodeAnum = randint(0,len(nodes))
                candnum  = range(len(nodes))
                candnum  = candnum[:nodeAnum] + candnum[nodeAnum+1:]
                tempidx  = randint(0,len(nodes)-1)
                nodeBnum = candnum[tempidx]
                
            def swap_nodes(nodeAnum, nodeBnum):
                def findNodes(root, nodeNum, nodeA=False, nodeB=False):
                    node = root
                    if nodeNum == nodeAnum:
                        nodeA = node
                    if nodeNum == nodeBnum:
                        nodeB = node
                    for i, child in enumerate(root['children']):
                        nodeNum = nodeNum + 1
                        (nodeA, nodeB, nodeNum) = findNodes(child, nodeNum, nodeA, nodeB)
                    return (nodeA, nodeB, nodeNum)
                
                (nodeA,nodeB,nodeNum) = findNodes(self.root, nodeNum=0)
                
                paramsA = nodeA['node'].params
                branchA = nodeA['node'].branch_length
                dataA = set(nodeA['node'].data)
                mainA = nodeA['main']
            
                nodeA['node'].params = nodeB['node'].params
                
                for dataid in list(dataA):
                    nodeA['node'].remove_datum(dataid)
                for dataid in nodeB['node'].data:
                    nodeA['node'].add_datum(dataid)
                    self.assignments[dataid] = nodeA['node']
                nodeA['main']=nodeB['main']
                
                nodeB['node'].params = paramsA
                dataB = set(nodeB['node'].data)
                
                for dataid in dataB:
                    nodeB['node'].remove_datum(dataid)
                for dataid in dataA:
                    nodeB['node'].add_datum(dataid)
                    self.assignments[dataid] = nodeB['node']
                nodeB['main']=mainA

            if empty_root:
                print 'checking alternative root'
                nodenum = []
                for ii, nn in enumerate(nodes):
                    if len(nodes[ii].data) > 0:
                        nodenum.append(ii)
                post_temp = zeros(len(nodenum))
                for idx, nodeBnum in enumerate(nodenum):
                    print 'nodeBnum', nodeBnum
                    print 'nodeAnum', nodeAnum
                    swap_nodes(nodeAnum, nodeBnum)
                    self.resample_sticks()
                    post_new = self.unnormalized_postertior()
                    post_temp[idx] = post_new
                    
                    if (post_new < post):
                        swap_nodes(nodeAnum,nodeBnum)
                        self.resample_sticks()
                        if nodeBnum == len(nodes)-1:
                            print 'forced swapping'
                            nodeBnum = post_temp.argmax() + 1
                            swap_nodes(nodeAnum,nodeBnum)
                            self.resample_sticks()
                            self.resample_node_params()
                            self.resample_stick_orders()
                    else:
                        print "successful swap!!!"
                        self.resample_node_params()
                        self.resample_stick_orders()
                        break
            else:   
                swap_nodes(nodeAnum,nodeBnum)
                self.resample_sticks()
                post_new = self.unnormalized_postertior()
      
                if (post_new < post):
                    swap_nodes(nodeAnum,nodeBnum)
                    self.resample_sticks()
                else:
                    print "successful swap!!!"
                    self.resample_node_params()
                    self.resample_stick_orders()
                                                

    def print_graph_full_logistic_different_branch_length(self, fh, base_width=5, min_width=300):
        edges   = sticks_to_edges(self.root['sticks'])
        weights = diff(hstack([0.0, edges]))
        if len(weights) > 0:
            root_mass = weights[0] * self.root['main']
        else: #in case there is a problem with weights, as observed in some runs
            root_mass = -1
        print >>fh, """graph: { title:            "TSSB Graph"  \
                                portsharing:      no            \
                                smanhattanedges:  yes           \
                                splines:          yes           \
                                equalydist:       yes           \
                                layout_algorithm: tree          \
                                node.fontname:    "helvR8"      \
                                node.height:      60            \
                                yspace:            20           \
                                xspace:            5 """
        print >>fh, """node: { label:"%d ~ Genotype """ %(len(self.root['node'].get_data())), \
            " ".join(map(lambda x: "%0.2f" %x, sigmoid(self.root['node'].params[range(8)]))), \
            "\n Branch Length: %0.2f Node Mass: %0.2f" % (self.root['node'].branch_length, root_mass), \
            """" title:"%s" width:%d}""" \
          %("X", min_width)
        def descend(root, name, mass):
            total   = 0.0
            edges   = sticks_to_edges(root['sticks'])
            weights = diff(hstack([0.0, edges]))
            for i, child in enumerate(root['children']):
                child_name = "%s-%d" % (name, i)
                child_mass = mass * weights[i] * child['main']
                print >>fh, """node: {  label:"%d ~ Genotype """ % (len(child['node'].get_data())), \
                " ".join(map(lambda x: "%0.2f" %x, sigmoid(child['node'].params[range(8)]))), \
                "\n Branch Length: %0.2f Node Mass: %0.2f" % (child['node'].branch_length, child_mass),\
                """" title:"%s" width:%d}""" \
                    %(child_name, min_width)
                print >>fh, """edge: { source:"%s" target:"%s" anchor:1}""" % (name, child_name)
                total += child_mass + descend(child, child_name, mass*weights[i] * (1.0 - child['main']))
            return total
        descend(self.root, 'X', 1)
        print >>fh, """}"""


    def tssb2igraph(self):
        
        edges   = sticks_to_edges(self.root['sticks'])
        weights = diff(hstack([0.0, edges]))
        if len(weights) > 0:
            root_mass = weights[0] * self.root['main']
        else: #in case there is a problem with weights, as observed in some runs
            root_mass = -1

        g = Graph(directed = True)
        g.add_vertex(name = "X", 
                     params = " ".join(map(lambda x: "%.15f" %x, 
                                           self.root['node'].params)), 
                     size = len(self.root['node'].get_data()), 
                     mass = root_mass, 
                     branch = self.root['node'].branch_length, 
                     members = " ".join(map(lambda x: "%s" %x, 
                                            self.root['node'].data)) ) 
        
        def descend(root, name, mass, g):
            total   = 0.0
            edges   = sticks_to_edges(root['sticks'])
            weights = diff(hstack([0.0, edges]))
            for i, child in enumerate(root['children']):
                child_name = "%s-%d" % (name, i)
                child_mass = mass * weights[i] * child['main']
                g.add_vertex(name = child_name, 
                             params = " ".join(map(lambda x: "%.15f" %x, 
                                                   child['node'].params)), 
                             size = len(child['node'].get_data()), 
                             mass = child_mass, 
                             branch = child['node'].branch_length,
                             members = " ".join(map(lambda x: "%s" %x, 
                                                    child['node'].data)) )
                g.add_edge(name, child_name, weight = child['node'].branch_length, 
                           Value = len(child['node'].get_data()))
                (tmp, g) = descend(child, 
                                   child_name, 
                                   mass*weights[i] * (1.0 - child['main']),
                                   g)
                
                total += child_mass + tmp
            return (total, g)

        (total,g) = descend(self.root, 'X', 1, g)
        return g

    def remove_empty_nodes(self):

        def descend(root):

            if len(root['children']) == 0:
                return

            while True:
                
                empty_nodes = filter(lambda i:
                                  len(root['children'][i]['node'].data) == 0, 
                                  range(len(root['children'])))

                if len(empty_nodes)==0:
                    break

                index = empty_nodes[0]
                
                cache_children = root['children'][index]['children']

                #root['children'][index]['node'].kill()

                del root['children'][index]

                if len(cache_children) == 0:
                    continue
                else:

                    temp1 = root['children'][:index]

                    temp2 = root['children'][index:]

                    root['children'] = temp1 + cache_children + temp2
                    root['sticks']   = zeros((len(root['children']),1))
                    
            for child in root['children']:
                descend(child)
            
                
        descend(self.root)
        self.resample_sticks()
    
    def deepest_node_depth(self):
        def descend(root, depth):
            if len(root['children']) == 0:
                return depth
            deepest = depth
            for i, child in enumerate(root['children']):
                hdepth = descend(child, depth+1)
                if deepest < hdepth:
                    deepest = hdepth
            return deepest 
        return descend(self.root, 1)
    
    def get_width_distribution(self):
        def descend(root, depth, width_vec):
            width_vec[depth-1] = width_vec[depth-1] + 1
            if len(root['children']) == 0:
                return width_vec
            for i, child in enumerate(root['children']):
                descend(child, depth+1, width_vec)
            return width_vec
        width_vec = zeros(self.max_depth)
        return descend(self.root, 1, width_vec)
    
    def get_weight_distribtuion(self):
        def descend(root, mass, depth, mass_vec):
            edges   = sticks_to_edges(root['sticks'])
            weights = diff(hstack([0.0, edges]))
            for i, child in enumerate(root['children']):
                mass_vec[depth] = mass_vec[depth] + mass * weights[i] * child['main']
                mass_vec = descend(child, mass*(1.0-child['main'])*weights[i],depth+1,mass_vec)
            return mass_vec 
        mass_vec = zeros(self.max_depth)
        edges   = sticks_to_edges(self.root['sticks'])
        weights = diff(hstack([0.0, edges]))
        if len(weights) > 0 :
            mass_vec[0] = weights[0] * self.root['main']
            return descend(self.root, 1.0-mass_vec[0], 1, mass_vec)
        else:
            return mass_vec
    
        
