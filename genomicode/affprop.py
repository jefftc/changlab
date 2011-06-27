"""

Functions:
cluster                 Run the affinity propagation algorithm.

set_preference_median   Set the preferences to be the median similarity.

"""
def cluster(similarity_matrix, dampfact=0.5, update_fn=None):
    # similarity_matrix is a NxN matrix of similarity scores.  More
    # similar points should have higher scores.  The score can be a
    # really small value if the two points are not connected.  Returns
    # a list of length N indicating the exemplar for each item.
    import itertools
    
    assert dampfact >= 0 and dampfact <= 1.0
    
    # Do some checking on the input matrix.
    nr = len(similarity_matrix)
    nc = len(similarity_matrix[0])
    assert nr == nc
    for x in similarity_matrix:
        assert len(x) == nc
    N = nr

    ## if not topology_matrix:
    ##     # Start with a fully connected graph.
    ##     topology_matrix = [[1]*N for i in range(N)]

    ## # Make sure topology matrix has the right dimensions.
    ## assert N == len(topology_matrix)
    ## assert N == len(topology_matrix[0])
    ## for x in topology_matrix:
    ##     assert len(x) == N

    ## # Make sure topology matrix is symmetric.
    ## for (i, j) in itertools.product(range(N), range(i, N)):
    ##     assert topology_matrix[i][j] == topology_matrix[j][i]

    ## topology = {}  # node -> list of connecting nodes.
    ## for i in range(N):
    ##     nodes = []
    ##     for j in range(N):
    ##         if topology_matrix[i][j]:
    ##             nodes.append(j)
    ##     topology[i] = nodes

    S = similarity_matrix
    # Availability matrix.
    A = [[0]*N for i in range(N)]
    # Responsibility matrix.
    R = [[0]*N for i in range(N)]

    df = dampfact
    num_iter = 0
    same_exemplars = 0
    exemplars = [None] * N
    #MAXITS = 2000
    #CONVITS = 200
    MAXITS = 100
    CONVITS = 5
    while num_iter < MAXITS and same_exemplars < CONVITS:
        num_iter += 1
        
        # Update the responsibility matrix.
        for i in range(N):
            # Find the largest and 2nd largest score.
            score1 = score_k = score2 = None
            for k in range(N):
                s = A[i][k] + S[i][k]
                if score1 is None or s > score1:
                    score1, score_k, score2 = s, k, score1
                elif score2 is None or s > score2:
                    score2 = s

            for k in range(N):
                if k != score_k:
                    R[i][k] = df*R[i][k] + (1-df)*(S[i][k]-score1)
                else:
                    R[i][k] = df*R[i][k] + (1-df)*(S[i][k]-score2)

        # Update the availability matrix.
        for k in range(N):
            scores = [max(0, R[i][k]) for i in range(N)]
            sum_scores = sum(scores)
            for i in range(N):
                s = sum_scores-scores[i]
                A[i][k] = df*A[i][k] + (1-df)*min(0, R[k][k]+s)
            s = sum_scores - scores[k]
            A[k][k] = df*A[k][k] + (1-df)*s

        # Calculate the exemplars.
        old_exemplars = exemplars
        exemplars = [None] * N
        for i in range(N):
            k_max = k_value = None
            for k in range(N):
                value = A[i][k] + R[i][k]
                if k_max is None or value > k_value:
                    k_max, k_value = k, value
            exemplars[i] = k_max
            
        # Count the number of changes in exemplars.
        changed = [int(x1 != x2) for (x1, x2) in zip(exemplars, old_exemplars)]
        num_changes = sum(changed)
        if num_changes == 0:
            same_exemplars += 1
        else:
            same_exemplars = 0
        if update_fn:
            update_fn(num_iter, num_changes, exemplars, A, R)

    clusters = _find_clusters(exemplars)
    return clusters, exemplars, S, A, R, num_iter, num_changes

def set_preference_median(similarity_matrix):
    # set s[k][k] to the median similarity.
    import jmath

    S = [x[:] for x in similarity_matrix]  # Make a copy.

    # Make sure S is a square matrix.
    N = len(S)
    for x in S:
        assert len(x) == N
    
    for i in range(len(S)):
        x = S[i][:]
        x.pop(i)
        S[i][i] = jmath.median(x)
    return S

def _find_clusters(exemplars):
    # Return a list of cluster IDs, where each cluster ID is an
    # integer from [0, N).  N is the number of clusters.
    
    # Figure out the topology of the graph.  Topology does not contain
    # nodes that point to themselves.
    topology = {}  # node -> list of connecting nodes.
    for i, j in enumerate(exemplars):
        if i == j:
            continue
        if i not in topology:
            topology[i] = []
        if j not in topology:
            topology[j] = []
        topology[i].append(j)
        topology[j].append(i)

    # Find the centers of each cluster.  A node is a cluster center
    # if:
    # 1.  It is connected to at least two other nodes.
    # 2.  Or it is a pair of nodes that are only connected to each
    #     other.  In this case, arbitrarily pick one to be the center.
    centers = []
    for node, partners in topology.iteritems():
        assert len(partners) > 0
        if len(partners) >= 2:
            centers.append(node)
        elif len(topology[partners[0]]) == 1 and partners[0] not in centers:
            centers.append(node)
    centers = sorted(centers)

    # Assign each node to a cluster, as defined by the centers.
    # 1.  If a node is a center, then set it to its own cluster.
    # 2.  If the node is connected to a center, then set it to that
    #     cluster.
    # 3.  Otherwise, set it to cluster 0.
    node2cluster = {}
    for node in topology:
        cluster = 0
        if node in centers:
            cluster = centers.index(node) + 1
        elif node in topology:
            for p in topology[node]:
                if p not in centers:
                    continue
                cluster = centers.index(p) + 1
                break
        node2cluster[node] = cluster

    clusters = [node2cluster.get(i, 0) for i in range(len(exemplars))]
    return clusters
    
        
def test_cluster():
    # 0,1 similar, 2,3,4 similar
    similarity_matrix = [
        [-3, -2, -4, -5, -6],
        [-3, -3, -6, -5, -4],
        [-4, -7, -3, -8, -3],
        [-2, -4, -2, -2, -1],
        [-5, -7, -5, -2, -2],
        ]
    x = cluster(similarity_matrix)
    print x

if __name__ == '__main__':
    test_cluster()
