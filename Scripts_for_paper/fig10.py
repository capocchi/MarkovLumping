

from multiprocessing import freeze_support
import platform
import pykov, sys, time

import numpy as np
import networkx as nx
from scipy.linalg import block_diag

from main_loop import *
from queuing import *

PLOT = True

def calculate_time(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"Function {func.__name__} took {(end_time - start_time):.5f} seconds to run.")
        return result
    return wrapper

def markov_chain_indexing_method(transition_matrix):
    """
    Applies Markov Chain Indexing Method to reduce the size of a transition matrix.
    
    Args:
        transition_matrix: np.array, the transition matrix of the Markov chain.
        num_indices: int, the number of indices to use for reduction.
        
    Returns:
        reduced_transition_matrix: np.array, the reduced transition matrix.
    """
    num_indices = len(transition_matrix)-1

    num_states = transition_matrix.shape[0]
    num_reduced_states = num_states // num_indices
    reduced_transition_matrix = np.zeros((num_reduced_states, num_reduced_states))
    
    for i in range(num_reduced_states):
        for j in range(num_reduced_states):
            for k in range(num_indices):
                for l in range(num_indices):
                    reduced_transition_matrix[i, j] += transition_matrix[i*num_indices + k, j*num_indices + l]
                    
    return reduced_transition_matrix

@calculate_time
def gth_reduction(P, tol=1e-10):
    """
    Implements the GTH reduction algorithm to reduce a Markov chain by eliminating
    transient states.
    
    Parameters:
        - P: numpy array representing the transition probability matrix of the Markov chain
        - tol: tolerance for convergence (default: 1e-10)
        
    Returns:
        - P_reduced: numpy array representing the reduced transition probability matrix
        - states: list of indices of the non-transient states in the reduced Markov chain
    """
    
    n = P.shape[0]
    D = np.diag(P.sum(axis=1))  # degree matrix
    
    # Compute initial GTH matrix
    M = D @ P
    M[np.diag_indices(n)] -= 1
    
    # Power method to find stationary distribution
    v = np.ones(n)
    v /= v.sum()
    while True:
        v_new = M @ v
        v_new /= v_new.sum()
        if np.max(np.abs(v_new - v)) < tol:
            break
        v = v_new
        
    # Identify non-transient states
    states = np.where(v > 0)[0]
    
    # Compute reduced transition probability matrix
    P_reduced = np.zeros((len(states), len(states)))
    for i, s in enumerate(states):
        for j, t in enumerate(states):
            P_reduced[i, j] = P[s, t] / v[t]
    
    return P_reduced, states.tolist()

@calculate_time
def reduce_graph(P, threshold=1e-5):
    """
    Réduit le graphe de la chaîne de Markov représentée par la matrice de transition P
    en supprimant les états dont l'importance est inférieure au seuil donné par threshold.
    """
    n = P.shape[0]  # nombre d'états de la chaîne de Markov
    G = nx.DiGraph(P)
    pr = nx.pagerank(G)  # calcule l'importance de chaque état
    to_remove = [i for i in range(n) if pr[i] < threshold]  # trouve les états à supprimer
    P_red = np.delete(np.delete(P, to_remove, axis=0), to_remove, axis=1)  # supprime les états
    return P_red

@calculate_time
def perron_frobenius_reduction(P, tol=1e-5):
    """
    Cette méthode permet de trouver une matrice de transition réduite pour une chaîne de Markov donnée.
    """
    n = P.shape[0]
    d = np.ones(n) / n
    while True:
        d_new = np.dot(d, P)
        if np.linalg.norm(d_new - d) < tol:
            break
        d = d_new
    Q = np.zeros((n-1, n-1))
    for i in range(n-1):
        for j in range(n-1):
            Q[i,j] = P[i,j] - d[i]*P[n-1,j] - P[i,n-1]*d[j] + d[n-1]*d[j]*P[n-1,n-1]
    return Q

@calculate_time
def perron_frobenius_reduction2(matrix):
    # Vérification que la matrice est carrée
    assert matrix.shape[0] == matrix.shape[1], "La matrice doit être carrée"
    
    # Vérification que la matrice est stochastique
    assert np.allclose(np.sum(matrix, axis=1), 1), "La matrice doit être stochastique"
    
    # Calcul des valeurs propres et vecteurs propres de la matrice
    eig_values, eig_vectors = np.linalg.eig(matrix.T)
    
    # Recherche de la valeur propre dominante
    max_eig_value_index = np.argmax(np.abs(eig_values))
    max_eig_value = eig_values[max_eig_value_index]
    max_eig_vector = eig_vectors[:, max_eig_value_index]
    
    # Normalisation du vecteur propre dominant
    max_eig_vector = max_eig_vector / np.sum(max_eig_vector)
    
    # Calcul de la matrice réduite ergodique
    reduced_matrix = np.outer(max_eig_vector, max_eig_vector)
    
    return reduced_matrix

@calculate_time
def qbd_reduce(P, Q):
    """
    Performs QBD reduction on a Markov chain with transition rate matrix Q,
    given the corresponding transition probability matrix P.

    Returns the reduced transition probability matrix and the group sizes.

    Parameters:
    - P: numpy array, the transition probability matrix
    - Q: numpy array, the transition rate matrix

    Returns:
    - P_reduced: numpy array, the reduced transition probability matrix
    - group_sizes: list of integers, the sizes of the groups of states
    """
    n = Q.shape[0]  # number of states
    R = np.diag(np.sum(Q, axis=1)) - Q  # the R matrix
    I = np.identity(n)

    # initialize the reduced matrices
    P_reduced = np.zeros((n, n))
    group_sizes = []

    # perform QBD reduction
    while True:
        # compute the quasi-diagonal blocks
        B = [I]
        for i in range(n):
            for j in range(n):
                if i != j and Q[i, j] != 0:
                    block = Q[i, j] * np.linalg.inv(R[i, i] * I - Q)
                    B.append(block)

        # compute the diagonal blocks and the group sizes
        D = []
        group_sizes = []
        i = 0
        while i < len(B):
            block = B[i]
            size = 1
            while i + size < len(B) and np.allclose(B[i + size], block):
                size += 1
            D.append(block)
            group_sizes.append(size)
            i += size

        # combine the blocks and check for convergence
        P_reduced = block_diag(*D)
        if np.allclose(P_reduced, P):
            break

        # compute the new Q and R matrices
        Q = np.sum([B[i] * R * B[i].T for i in range(len(B))], axis=0)
        R = np.diag(np.sum(Q, axis=1)) - Q

    return P_reduced, group_sizes

@calculate_time
def gerschgorin_reduction(P):
    # Compute the Gerschgorin disks
    n = P.shape[0]
    disks = [np.array([P[i,i]-np.abs(P[i,:]).sum(), P[i,i]+np.abs(P[i,:]).sum()]) for i in range(n)]
    
    # Group together disks that intersect
    groups = []
    for i in range(n):
        if i not in [group[-1] for group in groups]:
            groups.append([i])
        for j in range(i+1, n):
            if np.linalg.norm(disks[i] - disks[j]) <= disks[i][1] + disks[j][1]:
                if j not in [group[-1] for group in groups]:
                    groups[-1].append(j)
    
    # Construct the reduced transition probability matrix
    P_reduced = np.zeros((len(groups), len(groups)))
    for i, group_i in enumerate(groups):
        for j, group_j in enumerate(groups):
            for state_i in group_i:
                for state_j in group_j:
                    P_reduced[i,j] += P[state_i,state_j]
    
    # Normalize the rows of the reduced transition probability matrix
    row_sums = P_reduced.sum(axis=1)
    P_reduced = P_reduced / row_sums[:,np.newaxis]
    
    return P_reduced

@calculate_time
def pagerank_reduction(P, alpha=0.85, tol=1e-6, max_iter=10):
    """
    Reduce a Markov chain using the PageRank reduction algorithm.

    Parameters:
    -----------
    P : array-like, shape=(n_states, n_states)
        Transition probability matrix of the original Markov chain.
    alpha : float, optional
        Damping factor (default=0.85).
    tol : float, optional
        Tolerance for convergence (default=1e-6).
    max_iter : int, optional
        Maximum number of iterations (default=1000).

    Returns:
    --------
    R : array-like, shape=(n_reduced, n_reduced)
        Transition probability matrix of the reduced Markov chain.
    """
    n_states = P.shape[0]
    d = np.sum(P, axis=1)
    D_inv = np.diag(1 / d)
    Q = np.dot(D_inv, P)
    G = alpha * Q + (1 - alpha) * np.ones((n_states, n_states)) / n_states
    r = np.ones(n_states) / n_states
    for i in range(max_iter):
        r_new = np.dot(G, r)
        if np.linalg.norm(r_new - r, ord=1) < tol:
            break
        r = r_new
    mask = r > tol
    R = Q[mask][:, mask]
    return R

@calculate_time
def partial_sum_reduction(P, tol=1e-5):
    """
    Cette méthode utilise une approximation de la série géométrique pour calculer une matrice de transition réduite.
    """
    n = P.shape[0]
    Q = np.zeros((n-1, n-1))
    for i in range(n-1):
        for j in range(n-1):
            if i == j:
                Q[i,j] = 1 - P[n-1,n-1]
            else:
                Q[i,j] = P[i,j] / (1 - P[n-1,n-1])
    return Q

def mm1_queue_stats(lambd, mu, c):
    """Calcule les statistiques de la queue M/M/1
    
    Args:
        lambd (float): le taux d'arrivée moyen des clients
        mu (float): le taux de service moyen
        c (float): la capacité du système
    
    Returns:
        Tuple[float, float, float, float]: le temps d'attente moyen, le temps de service moyen, 
        le taux de service et le taux d'abandon
    """
    rho = lambd / mu  # facteur de charge
    p0 = 1 - rho  # probabilité d'état 0
    Lq = rho**2 / (1 - rho)  # nombre moyen de clients dans la file d'attente
    Wq = Lq / lambd  # temps moyen d'attente dans la file d'attente
    W = Wq + 1 / mu  # temps moyen de service (attente +

def generate_transition_matrix(arrival_rate, service_rate, num_servers):
    # Calculate the arrival rate per server
    arrival_rate_per_server = arrival_rate / num_servers
    
    # Calculate the transition rates for each state
    transition_rates = np.zeros((num_servers + 1, num_servers + 1))
    for i in range(num_servers):
        transition_rates[i][i+1] = arrival_rate_per_server
        transition_rates[i+1][i] = min(i, num_servers - i) * service_rate
        transition_rates[i][i] = -(arrival_rate_per_server + (i * service_rate))
    
    # Calculate the diagonal elements of the transition matrix
    diagonal = np.abs(transition_rates).sum(axis=1)
    
    # Set the diagonal elements of the transition matrix
    transition_matrix = np.diag(diagonal)
    
    # Set the off-diagonal elements of the transition matrix
    for i in range(num_servers + 1):
        for j in range(num_servers + 1):
            if i != j:
                transition_matrix[i][j] = transition_rates[i][j]
    
    # Normalize the rows of the transition matrix
    row_sums = transition_matrix.sum(axis=1)
    transition_matrix = transition_matrix / row_sums[:, np.newaxis]
    
    return transition_matrix

def generate_transition_matrix2(arrival_rate, service_rate, num_servers):
    rho = arrival_rate / (num_servers * service_rate)
    p = np.zeros((num_servers+1, num_servers+1))
    p[0,0] = 1 - rho
    for i in range(1, num_servers+1):
        p[i,i-1] = rho * p[i-1,i-1] / i
        p[i,i] = 1 - sum(p[i,:-1])
        p[i,i+1:] = rho * p[i,i:-1] / (i+1)
    return p

def generate_complex_transition_matrix(n, k):
    # Créer une matrice n x n avec des valeurs aléatoires entre 0 et 1
    M = np.random.rand(n, n)
    
    # Définir les k entrées les plus grandes pour chaque ligne à 1, les autres à 0
    for i in range(n):
        top_k_indices = M[i, :].argsort()[-k:]
        M[i, :] = 0
        M[i, top_k_indices] = 1
    
    # Normaliser chaque ligne pour avoir une somme de 1
    for i in range(n):
        row_sum = M[i, :].sum()
        M[i, :] = M[i, :] / row_sum
        
    return M

def getStat3(queue_transition_matrix, arrival_rate, service_rate, num_servers):
    
    # Compute steady-state probabilities
    pi = np.linalg.solve(np.transpose(queue_transition_matrix) - np.eye(num_servers+1), np.zeros(num_servers+1) + 1)
    pi /= sum(pi)

    # Compute average waiting time, service time, service rate, and abandon rate
    avg_waiting_time = sum([i * pi[i] for i in range(1, num_servers+1)]) / (num_servers * service_rate - arrival_rate)
    avg_service_time = 1 / service_rate
    service_rate_effective = num_servers * service_rate * (1 - pi[0])
    abandon_rate = arrival_rate * pi[num_servers] / service_rate_effective

    print("Average waiting time: {:.2f}".format(avg_waiting_time))
    print("Average service time: {:.2f}".format(avg_service_time))
    print("Service rate: {:.2f}".format(service_rate_effective))
    print("Abandon rate: {:.2f}".format(abandon_rate))

    return (avg_waiting_time, avg_service_time, service_rate_effective, abandon_rate)

def plot1(X, normalized_QUEUE_WAITING_TIME_ORIGINAL,normalized_QUEUE_WAITING_TIME_BESTA,
                normalized_QUEUE_WAITING_TIME_M2,normalized_QUEUE_WAITING_TIME_M3,
                normalized_QUEUE_WAITING_TIME_M4,normalized_kl):
    
     # # Créer une figure avec deux sous-graphiques
    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(nrows=5, ncols=1)
    
    #ax1.plot(X,normalized_kl[0], label="KL", marker="o")
    ax1.plot(X,normalized_QUEUE_WAITING_TIME_ORIGINAL, label="Original", marker=".")
    ax1.plot(X,normalized_QUEUE_WAITING_TIME_BESTA, label="Besta", marker="1")
    #ax1.plot(X,normalized_QUEUE_WAITING_TIME_M1, label="perron_frobenius", marker="*")
    ax1.plot(X,normalized_QUEUE_WAITING_TIME_M2, label="gerschgorin", marker="o")
    ax1.plot(X,normalized_QUEUE_WAITING_TIME_M3, label="partial_sum", marker=".")
    ax1.plot(X,normalized_QUEUE_WAITING_TIME_M4, label="page_rank", marker="*")
    ax1.set_title("Average Waiting time")
    ax1.tick_params(rotation=45)
    # show a legend on the plot
    #ax1.legend()
    ax1.grid()
    ax1.set_ylabel("Normalized Amplitude")

    #ax2.plot(X,normalized_kl[0], label="KL", marker="o")
    ax2.plot(X,normalized_QUEUE_SERVICE_TIME_ORIGINAL, marker=".")
    ax2.plot(X,normalized_QUEUE_SERVICE_TIME_BESTA, marker="1")
    #ax2.plot(X,normalized_QUEUE_SERVICE_TIME_M1[0], marker="*")
    ax2.plot(X,normalized_QUEUE_SERVICE_TIME_M2, marker="o")
    ax2.plot(X,normalized_QUEUE_SERVICE_TIME_M3, marker=".")
    ax2.plot(X,normalized_QUEUE_SERVICE_TIME_M4, marker="*")
    ax2.set_title("Average Service time")
    ax2.tick_params(rotation=45)
    # show a legend on the plot
    #ax2.legend()
    ax2.grid()
    ax2.set_ylabel("Normalized Amplitude")

    #ax3.plot(X,normalized_kl[0], label="KL", marker="o")
    ax3.plot(X,normalized_QUEUE_SERVICE_RATE_ORIGINAL, marker=".")
    ax3.plot(X,normalized_QUEUE_SERVICE_RATE_BESTA, marker="1")
    #ax3.plot(X,normalized_QUEUE_SERVICE_RATE_M1[0], marker="*")
    ax3.plot(X,normalized_QUEUE_SERVICE_RATE_M2,  marker="o")
    ax3.plot(X,normalized_QUEUE_SERVICE_RATE_M3, marker=".")
    ax3.plot(X,normalized_QUEUE_SERVICE_RATE_M4, marker="*")
    ax3.set_title("Service Rate")
    ax3.tick_params(rotation=45)
    ax3.grid()
    ax3.set_ylabel("Normalized Amplitude")

    #ax4.plot(X,normalized_kl[0], label="KL", marker="o")
    ax4.plot(X,normalized_QUEUE_ABANDON_RATE_ORIGINAL, marker=".")
    ax4.plot(X,normalized_QUEUE_ABANDON_RATE_BESTA, marker="1")
    #ax4.plot(X,normalized_QUEUE_ABANDON_RATE_M1[0], marker="*")
    ax4.plot(X,normalized_QUEUE_ABANDON_RATE_M2, marker="o")
    ax4.plot(X,normalized_QUEUE_ABANDON_RATE_M3, marker=".")
    ax4.plot(X,normalized_QUEUE_ABANDON_RATE_M4, marker="*")
    ax4.set_title("Abandon Rate")
    ax4.tick_params(rotation=45)
    ax4.grid()
    ax4.set_ylabel("Normalized Amplitude")

    ax5.plot(X,normalized_kl, label="KL", marker="o")
    ax5.set_title("KL")
    ax5.tick_params(rotation=45)
    ax5.grid()
    ax5.set_ylabel("Normalized Amplitude")

    # # Introduire un espace entre les subplots
    fig.subplots_adjust(hspace=0.5)

    fig.legend(loc='outside right')

    plt.plot()

    # function to show the plot
    plt.show()

def plot2(X,Y=[]):
    
    # Plot the data on each subplot
    # axes[0, 0].scatter(X, normalized_QUEUE_ABANDON_RATE_ORIGINAL[0], marker="o", label="Original")
    
    axes[0, 0].scatter(X, QUEUE_ABANDON_RATE_BESTA[0], marker="*", label="BESTA")
    axes[0, 0].scatter(X, QUEUE_ABANDON_RATE_M2[0], marker="h", label="Gerschgorin")
    axes[0, 0].scatter(X, QUEUE_ABANDON_RATE_M3[0], marker="D", label = "Partial Sum")
    axes[0, 0].scatter(X, QUEUE_ABANDON_RATE_M4[0], marker="p", label="Page Rank")
    axes[0, 0].set_title("Abandon Rate")
    axes[0, 0].tick_params(rotation=45)
    axes[0, 0].legend()

    # axes[0, 1].scatter(X, normalized_QUEUE_SERVICE_RATE_ORIGINAL[0], marker="o", label="Original")
    axes[0, 1].scatter(X, QUEUE_SERVICE_RATE_BESTA[0], marker="*", label="BESTA")
    axes[0, 1].scatter(X, QUEUE_SERVICE_RATE_M2[0], marker="h", label="Gerschgorin")
    axes[0, 1].scatter(X, QUEUE_SERVICE_RATE_M3[0], marker="D", label = "Partial Sum")
    axes[0, 1].scatter(X, QUEUE_SERVICE_RATE_M4[0], marker="p", label="Page Rank")
    axes[0, 1].set_title("Service Rate")
    axes[0, 1].tick_params(rotation=45)
    # axes[0, 1].legend()

    # axes[1, 0].scatter(X, normalized_QUEUE_SERVICE_TIME_ORIGINAL[0], marker="o", label="Original")
    axes[1, 0].scatter(X, QUEUE_SERVICE_TIME_BESTA[0], marker="*", label="BESTA")
    axes[1, 0].scatter(X, QUEUE_SERVICE_TIME_M2[0], marker="h", label="Gerschgorin")
    axes[1, 0].scatter(X, QUEUE_SERVICE_TIME_M3[0], marker="D", label = "Partial Sum")
    axes[1, 0].scatter(X, QUEUE_SERVICE_TIME_M4[0], marker="p", label="Page Rank")
    axes[1, 0].set_title("Service Time")
    axes[1, 0].tick_params(rotation=45)
    # axes[1, 0].legend()

    # axes[1, 1].scatter(X, normalized_QUEUE_WAITING_TIME_ORIGINAL[0], marker="o", label="Original")
    axes[1, 1].scatter(X, QUEUE_WAITING_TIME_BESTA[0], marker="*", label="BESTA")
    axes[1, 1].scatter(X, QUEUE_WAITING_TIME_M2[0], marker="h", label="Gerschgorin")
    axes[1, 1].scatter(X, QUEUE_WAITING_TIME_M3[0], marker="D", label = "Partial Sum")
    axes[1, 1].scatter(X, QUEUE_WAITING_TIME_M4[0], marker="p", label="Page Rank")
    axes[1, 1].set_title("Waiting Time")
    axes[1, 1].tick_params(rotation=45)
    # axes[1, 1].legend()

    plt.ylabel("Normalized Amplitude")
    #plt.vlines(range(0,10), ymin=-1, ymax=1, colors='gray', linestyles='dashed')

    # Adjust the spacing between subplots
    fig.subplots_adjust(hspace=0.3, wspace=0.3)

    # function to show the plot
    plt.show()

    return transition_matrix

if __name__ == '__main__':
    """ to use : python fig10.py 19 40 20
    """
    from Queue import Queue

    # for Windows support of tqdm
    if platform == "win32":
        freeze_support()

    ### size of queue
    try:
        n = sys.argv[1]
        sr = sys.argv[2]
        ar = sys.argv[3]
    except:
        sys.exit()
    else:
        
        num_servers = int(n)
        service_rate = float(sr)
        arrival_rate = float(ar)
    
        # Le temps moyen d'attente dans la file d'attente
        #mean_waiting_time = (1 / (service_rate - arrival_rate)) * (1 - (num_servers * (arrival_rate / service_rate))**num_servers * (1 - (arrival_rate / service_rate))**(num_servers + 1))

        # facteur de charge 
        #rho = arrival_rate / (num_servers * service_rate)

        #queue_transition_matrix = generate_mm1_queue_chain(arrival_rate,service_rate,num_servers)
        
        queue_transition_matrix= generate_queue_markov_chain(arrival_rate,service_rate,num_servers)

        # queue_transition_matrix=  generate_transition_matrix2(arrival_rate,service_rate,num_servers)

        methods = { #'perron_frobenius':perron_frobenius_reduction2, 
                    'gerschgorin':gerschgorin_reduction, 
                    #'GTH':gth_reduction,
                    'partial_sum':partial_sum_reduction,
                    'chain_indexing':markov_chain_indexing_method}
                    # 'page_rank':pagerank_reduction}

        ### get queue stat
        print("---------------------- original constants")
        avg_waiting_time,avg_service_time,service_rate,abandon_rate = getStat3(queue_transition_matrix, arrival_rate, service_rate, num_servers)
        QUEUE_WAITING_TIME_ORIGINAL = [avg_waiting_time]
        QUEUE_SERVICE_TIME_ORIGINAL = [avg_service_time]
        QUEUE_SERVICE_RATE_ORIGINAL = [service_rate]
        QUEUE_ABANDON_RATE_ORIGINAL = [abandon_rate]

        #exit()

        ### first lumping with BESTA
        P = pykov.Chain()
        for i,row in enumerate(queue_transition_matrix):
            for j,val in enumerate(row):
                P[(str(i),str(j))]=val
        
        S = tuple(sorted(list(P.states())))
       

        # starting time
        start1 = time.time()
        p,Q = next(getMFTPAnalysis3(S,P))
        # ending time
        end1 = time.time()
         # total time taken
        print(f"Function BESTA took {(end1 - start1):.5f} seconds to run.")
        
        n = len(S)
       
        A = chainToNPArray(Q)

        ### get queue stat
        print(f"---------------------------------------------------------- {n}x{n}")
        print("---------------------------------------- BESTA")
        avg_waiting_time,avg_service_time,service_rate,abandon_rate = getStat3(A, arrival_rate, service_rate, num_servers-1)
        QUEUE_WAITING_TIME_BESTA = [avg_waiting_time]
        QUEUE_SERVICE_TIME_BESTA = [avg_service_time]
        QUEUE_SERVICE_RATE_BESTA = [service_rate]
        QUEUE_ABANDON_RATE_BESTA = [abandon_rate]

        for k,v in methods.items():
            print(f"---------------------- {k}\n")
            awt,ast,sr,ar = getStat(v(A),True)
            
            if awt < 0: awt = 0
            if ast < 0: ast = 0
            if sr < 0: sr = 0
            if ar < 0 : ar = 0
            
            if k == "perron_frobenius":
                QUEUE_WAITING_TIME_M1 = [awt]
                QUEUE_SERVICE_TIME_M1 = [ast]
                QUEUE_SERVICE_RATE_M1 = [sr]
                QUEUE_ABANDON_RATE_M1 = [ar]
            elif k == "gerschgorin":
                QUEUE_WAITING_TIME_M2 = [awt]
                QUEUE_SERVICE_TIME_M2 = [ast]
                QUEUE_SERVICE_RATE_M2 = [sr]
                QUEUE_ABANDON_RATE_M2 = [ar]
            elif k == "partial_sum":
                QUEUE_WAITING_TIME_M3 = [awt]
                QUEUE_SERVICE_TIME_M3 = [ast]
                QUEUE_SERVICE_RATE_M3 = [sr]
                QUEUE_ABANDON_RATE_M3 = [ar]
            elif k == "page_rank":
                QUEUE_WAITING_TIME_M4 = [awt]
                QUEUE_SERVICE_TIME_M4 = [ast]
                QUEUE_SERVICE_RATE_M4 = [sr]
                QUEUE_ABANDON_RATE_M4 = [ar]
            elif k == "chain_indexing":
                QUEUE_WAITING_TIME_M5 = [awt]
                QUEUE_SERVICE_TIME_M5 = [ast]
                QUEUE_SERVICE_RATE_M5 = [sr]
                QUEUE_ABANDON_RATE_M5 = [ar]
            else:
                exit()
                        
        if PLOT:
            X = [n]

        ### condition for table 2
        cond = "n>2"

        ### condition for table3
        # kl=new_kl=diff=0.0
        # cond = "new_kl <= kl*(1+0.5) or kl==0.0"

        ### condition for table4
        # service_rate = 0.0
        # cond = "abs(service_rate)-QUEUE_SERVICE_RATE_ORIGINAL[0]) > abs(QUEUE_SERVICE_RATE_ORIGINAL[0])*(1+0.1) or service_rate==0.0"

        n = len(A)

        ### stopping condition
        while(eval(cond)):
            
            # starting time
            start2 = time.time()

            ### P must be ergotic i.e. the transition matrix must be irreducible and acyclic.            
            G = nx.DiGraph(list(P.keys()), directed=True)
            nx.strongly_connected_components(G)
            assert nx.is_strongly_connected(G) and nx.is_aperiodic(G), f"Matrix is not ergotic!"
            
            ### just for stdout
            d = trace(n,0.0,p)
            
            ### P transition matrix of the Markoc chain to lump
            P = pykov.Chain()
            for k,v in Q.items():
                P[(d[k[0]],d[k[1]])] = v

            A = chainToNPArray(P)
            ### get queue stat
            print("---------------------- BESTA")
            avg_waiting_time,avg_service_time,service_rate,abandon_rate = getStat(A, True)
            
            print(abs(service_rate), abs(QUEUE_SERVICE_RATE_ORIGINAL[0]),abs(QUEUE_SERVICE_RATE_ORIGINAL[0])*(1+0.1))

            QUEUE_WAITING_TIME_BESTA.append(avg_waiting_time)
            QUEUE_SERVICE_TIME_BESTA.append(avg_service_time)
            QUEUE_SERVICE_RATE_BESTA.append(service_rate)
            QUEUE_ABANDON_RATE_BESTA.append(abandon_rate)

            for k,v in methods.items():
                print(f"---------------------- {k}\n")
                awt,ast,sr,ar = getStat(v(A),True)
                
                if awt < 0: awt = 0
                if ast < 0: ast = 0
                if sr < 0: sr = 0
                if ar < 0 : ar = 0

                if k == "perron_frobenius":
                    QUEUE_WAITING_TIME_M1.append(awt)
                    QUEUE_SERVICE_TIME_M1.append(ast)
                    QUEUE_SERVICE_RATE_M1.append(sr)
                    QUEUE_ABANDON_RATE_M1.append(ar)
                elif k == "gerschgorin":
                    QUEUE_WAITING_TIME_M2.append(awt)
                    QUEUE_SERVICE_TIME_M2.append(ast)
                    QUEUE_SERVICE_RATE_M2.append(sr)
                    QUEUE_ABANDON_RATE_M2.append(ar)
                elif k == "partial_sum":
                    QUEUE_WAITING_TIME_M3.append(awt)
                    QUEUE_SERVICE_TIME_M3.append(ast)
                    QUEUE_SERVICE_RATE_M3.append(sr)
                    QUEUE_ABANDON_RATE_M3.append(ar)
                elif k == "page_rank":
                    QUEUE_WAITING_TIME_M4.append(awt)
                    QUEUE_SERVICE_TIME_M4.append(ast)
                    QUEUE_SERVICE_RATE_M4.append(sr)
                    QUEUE_ABANDON_RATE_M4.append(ar)
                elif k == "chain_indexing":
                    QUEUE_WAITING_TIME_M5.append(awt)
                    QUEUE_SERVICE_TIME_M5.append(ast)
                    QUEUE_SERVICE_RATE_M5.append(sr)
                    QUEUE_ABANDON_RATE_M5.append(ar)
                else:
                    exit()
                    
            ### set of states
            S = tuple(sorted(P.states()))

            ### mfpt analisys (REDUCED function that call FOP)
            p,Q = next(getMFTPAnalysis3(S,P))

            ### update variables
            n = len(S)

            if PLOT:
                X.append(n)
            

            if WRITE_FILE:
                fn = os.path.join(os.pardir,'Matrix',f"{n}x{n}.dat")
                if os.path.exists(fn):
                    os.remove(fn)
                f = open(fn,'w')
                for k,v in dict(P).items():
                    f.write(f"{k[0]} {k[1]} {v} \n")
                f.close()
         
        ### just for stdout
        # trace(n,kl,p)
        A = chainToNPArray(P)
      
        if PLOT:
            X = list(map(lambda a : f"{a} x {a}", map(str,X)))

            QUEUE_WAITING_TIME_ORIGINAL = QUEUE_WAITING_TIME_ORIGINAL*len(X)
            QUEUE_SERVICE_TIME_ORIGINAL = QUEUE_SERVICE_TIME_ORIGINAL*len(X)
            QUEUE_SERVICE_RATE_ORIGINAL = QUEUE_SERVICE_RATE_ORIGINAL*len(X)
            QUEUE_ABANDON_RATE_ORIGINAL = QUEUE_ABANDON_RATE_ORIGINAL*len(X)

            normalized_QUEUE_WAITING_TIME_ORIGINAL = list(map(lambda a: a/max(QUEUE_WAITING_TIME_ORIGINAL), QUEUE_WAITING_TIME_ORIGINAL)) #preprocessing.normalize([QUEUE_WAITING_TIME_ORIGINAL])
            normalized_QUEUE_WAITING_TIME_BESTA = list(map(lambda a: a/max(QUEUE_WAITING_TIME_BESTA), QUEUE_WAITING_TIME_BESTA))#preprocessing.normalize([QUEUE_WAITING_TIME_BESTA])
            normalized_QUEUE_SERVICE_TIME_ORIGINAL = list(map(lambda a: a/max(QUEUE_SERVICE_TIME_ORIGINAL), QUEUE_SERVICE_TIME_ORIGINAL))#preprocessing.normalize([QUEUE_SERVICE_TIME_ORIGINAL])
            normalized_QUEUE_SERVICE_TIME_BESTA = list(map(lambda a: a/max(QUEUE_SERVICE_TIME_BESTA), QUEUE_SERVICE_TIME_BESTA)) #preprocessing.normalize([QUEUE_SERVICE_TIME_BESTA])
            normalized_QUEUE_SERVICE_RATE_ORIGINAL = list(map(lambda a: a/max(QUEUE_SERVICE_RATE_ORIGINAL), QUEUE_SERVICE_RATE_ORIGINAL)) #preprocessing.normalize([QUEUE_SERVICE_RATE_ORIGINAL])
            normalized_QUEUE_SERVICE_RATE_BESTA = list(map(lambda a: a/max(QUEUE_SERVICE_RATE_BESTA), QUEUE_SERVICE_RATE_BESTA)) #preprocessing.normalize([QUEUE_SERVICE_RATE_BESTA])
            normalized_QUEUE_ABANDON_RATE_ORIGINAL = list(map(lambda a: a/max(QUEUE_ABANDON_RATE_ORIGINAL), QUEUE_ABANDON_RATE_ORIGINAL)) #preprocessing.normalize([QUEUE_ABANDON_RATE_ORIGINAL])
            normalized_QUEUE_ABANDON_RATE_BESTA = list(map(lambda a: a/max(QUEUE_ABANDON_RATE_BESTA), QUEUE_ABANDON_RATE_BESTA)) #preprocessing.normalize([QUEUE_ABANDON_RATE_BESTA])

            for k in methods:
                if k == "perron_frobenius":
                    normalized_QUEUE_WAITING_TIME_M1 = list(map(lambda a: a/max(QUEUE_WAITING_TIME_M1), QUEUE_WAITING_TIME_M1)) #preprocessing.normalize([QUEUE_WAITING_TIME_M1])
                    normalized_QUEUE_SERVICE_TIME_M1 = list(map(lambda a: a/max(QUEUE_SERVICE_TIME_M1), QUEUE_SERVICE_TIME_M1)) #preprocessing.normalize([QUEUE_SERVICE_TIME_M1])
                    normalized_QUEUE_SERVICE_RATE_M1 = list(map(lambda a: a/max(QUEUE_SERVICE_RATE_M1), QUEUE_SERVICE_RATE_M1)) #preprocessing.normalize([QUEUE_SERVICE_RATE_M1])
                    normalized_QUEUE_ABANDON_RATE_M1 = list(map(lambda a: a/max(QUEUE_ABANDON_RATE_M1), QUEUE_ABANDON_RATE_M1)) #preprocessing.normalize([QUEUE_ABANDON_RATE_M1])
                elif k == "gerschgorin":
                    normalized_QUEUE_WAITING_TIME_M2 = list(map(lambda a: a/max(QUEUE_WAITING_TIME_M2), QUEUE_WAITING_TIME_M2)) #preprocessing.normalize([QUEUE_WAITING_TIME_M2])
                    normalized_QUEUE_SERVICE_TIME_M2 = list(map(lambda a: a/max(QUEUE_SERVICE_TIME_M2), QUEUE_SERVICE_TIME_M2)) #preprocessing.normalize([QUEUE_SERVICE_TIME_M2])
                    normalized_QUEUE_SERVICE_RATE_M2 = list(map(lambda a: a/max(QUEUE_SERVICE_RATE_M2), QUEUE_SERVICE_RATE_M2)) #preprocessing.normalize([QUEUE_SERVICE_RATE_M2])
                    normalized_QUEUE_ABANDON_RATE_M2 = list(map(lambda a: a/max(QUEUE_ABANDON_RATE_M2), QUEUE_ABANDON_RATE_M2)) #preprocessing.normalize([QUEUE_ABANDON_RATE_M2])
                elif k == "partial_sum":
                    normalized_QUEUE_WAITING_TIME_M3 = list(map(lambda a: a/max(QUEUE_WAITING_TIME_M3), QUEUE_WAITING_TIME_M3))  #preprocessing.normalize([QUEUE_WAITING_TIME_M3])
                    normalized_QUEUE_SERVICE_TIME_M3 = list(map(lambda a: a/max(QUEUE_SERVICE_TIME_M3), QUEUE_SERVICE_TIME_M3))#preprocessing.normalize([QUEUE_SERVICE_TIME_M3])
                    normalized_QUEUE_SERVICE_RATE_M3 = list(map(lambda a: a/max(QUEUE_SERVICE_RATE_M3), QUEUE_SERVICE_RATE_M3)) #preprocessing.normalize([QUEUE_SERVICE_RATE_M3])
                    normalized_QUEUE_ABANDON_RATE_M3 = list(map(lambda a: a/max(QUEUE_ABANDON_RATE_M3), QUEUE_ABANDON_RATE_M3))#preprocessing.normalize([QUEUE_ABANDON_RATE_M3])
                elif k == "page_rank":
                    normalized_QUEUE_WAITING_TIME_M4 = list(map(lambda a: a/max(QUEUE_WAITING_TIME_M4), QUEUE_WAITING_TIME_M4)) #preprocessing.normalize([QUEUE_WAITING_TIME_M4])
                    normalized_QUEUE_SERVICE_TIME_M4 = list(map(lambda a: a/max(QUEUE_SERVICE_TIME_M4), QUEUE_SERVICE_TIME_M4))#preprocessing.normalize([QUEUE_SERVICE_TIME_M4])
                    normalized_QUEUE_SERVICE_RATE_M4 = list(map(lambda a: a/max(QUEUE_SERVICE_RATE_M4), QUEUE_SERVICE_RATE_M4)) #preprocessing.normalize([QUEUE_SERVICE_RATE_M4])
                    normalized_QUEUE_ABANDON_RATE_M4 = list(map(lambda a: a/max(QUEUE_ABANDON_RATE_M4), QUEUE_ABANDON_RATE_M4))#preprocessing.normalize([QUEUE_ABANDON_RATE_M4])
                elif k == "chain_indexing":
                    normalized_QUEUE_WAITING_TIME_M5 = list(map(lambda a: a/max(QUEUE_WAITING_TIME_M5), QUEUE_WAITING_TIME_M5)) #preprocessing.normalize([QUEUE_WAITING_TIME_M5])
                    normalized_QUEUE_SERVICE_TIME_M5 = list(map(lambda a: a/max(QUEUE_SERVICE_TIME_M5), QUEUE_SERVICE_TIME_M5))#preprocessing.normalize([QUEUE_SERVICE_TIME_M5])
                    normalized_QUEUE_SERVICE_RATE_M5 = list(map(lambda a: a/max(QUEUE_SERVICE_RATE_M5), QUEUE_SERVICE_RATE_M5)) #preprocessing.normalize([QUEUE_SERVICE_RATE_M5])
                    normalized_QUEUE_ABANDON_RATE_M5 = list(map(lambda a: a/max(QUEUE_ABANDON_RATE_M5), QUEUE_ABANDON_RATE_M5))#preprocessing.normalize([QUEUE_ABANDON_RATE_M5])
                else:
                    pass

            fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))

            # Plot the data on each subplot
            axes[0, 0].plot(X, normalized_QUEUE_ABANDON_RATE_ORIGINAL, label="Original")
            axes[0, 0].scatter(X, normalized_QUEUE_ABANDON_RATE_BESTA, marker="*", label = "BESTA")
            axes[0, 0].scatter(X, normalized_QUEUE_ABANDON_RATE_M2, marker="h", label = "Gerschgorin")
            axes[0, 0].scatter(X, normalized_QUEUE_ABANDON_RATE_M3, marker="D", label = "Partial Sum")
            axes[0, 0].scatter(X, normalized_QUEUE_ABANDON_RATE_M5, marker="p", label = "Chain Indexing")
            axes[0, 0].set_title("Abandon Rate")
            axes[0, 0].tick_params(rotation=45)
            axes[0, 0].legend()
            axes[0, 0].set_ylabel("Normalized Amplitude")
            axes[0, 0].grid(axis='x', color='#D3D3D3')

            axes[0, 1].plot(X, normalized_QUEUE_SERVICE_RATE_ORIGINAL, label="Original")
            axes[0, 1].scatter(X, normalized_QUEUE_SERVICE_RATE_BESTA, marker="*", label = "BESTA")
            axes[0, 1].scatter(X, normalized_QUEUE_SERVICE_RATE_M2, marker="h", label = "Gerschgorin")
            axes[0, 1].scatter(X, normalized_QUEUE_SERVICE_RATE_M3, marker="D", label = "Partial Sum")
            axes[0, 1].scatter(X, normalized_QUEUE_SERVICE_RATE_M5, marker="p", label = "Chain Indexing")
            axes[0, 1].set_title("Service Rate")
            axes[0, 1].tick_params(rotation=45)
            axes[0, 1].set_ylabel("Normalized Amplitude")
            axes[0, 1].grid(axis='x', color='#D3D3D3')

            axes[1, 0].plot(X, normalized_QUEUE_SERVICE_TIME_ORIGINAL, label="Original")
            axes[1, 0].scatter(X, normalized_QUEUE_SERVICE_TIME_BESTA, marker="*", label = "BESTA")
            axes[1, 0].scatter(X, normalized_QUEUE_SERVICE_TIME_M2, marker="h", label = "Gerschgorin")
            axes[1, 0].scatter(X, normalized_QUEUE_SERVICE_TIME_M3, marker="D", label = "Partial Sum")
            axes[1, 0].scatter(X, normalized_QUEUE_SERVICE_TIME_M5, marker="p", label = "Chain Indexing")
            axes[1, 0].set_title("Average Service Time")
            axes[1, 0].tick_params(rotation=45)
            axes[1, 0].set_ylabel("Normalized Amplitude")
            axes[1, 0].grid(axis='x', color='#D3D3D3')

            axes[1, 1].plot(X, normalized_QUEUE_WAITING_TIME_ORIGINAL, label="Original")
            axes[1, 1].scatter(X, normalized_QUEUE_WAITING_TIME_BESTA, marker="*", label = "BESTA")
            axes[1, 1].scatter(X, normalized_QUEUE_WAITING_TIME_M2, marker="h", label = "Gerschgorin")
            axes[1, 1].scatter(X, normalized_QUEUE_WAITING_TIME_M3, marker="D", label = "Partial Sum")
            axes[1, 1].scatter(X, normalized_QUEUE_WAITING_TIME_M5, marker="p", label ="Chain Indexing")
            axes[1, 1].set_title("Average Waiting Time")
            axes[1, 1].tick_params(rotation=45)
            axes[1, 1].set_ylabel("Normalized Amplitude")
            axes[1, 1].grid(axis='x', color='#D3D3D3')

            plt.ylabel("Normalized Amplitude")
            #plt.vlines(range(0,10), ymin=-1, ymax=1, colors='gray', linestyles='dashed')

            # Adjust the spacing between subplots
            fig.subplots_adjust(hspace=0.3, wspace=0.3)

            plt.savefig("compa.png", dpi=600)

            # function to show the plot
            plt.show()
        
          
        from sklearn.metrics import mean_squared_error
        from math import sqrt
        from collections import OrderedDict

        original_abandon = QUEUE_ABANDON_RATE_ORIGINAL[0]
        original_sr = QUEUE_SERVICE_RATE_ORIGINAL[0]
        original_st = QUEUE_SERVICE_TIME_ORIGINAL[0]
        original_wt = QUEUE_WAITING_TIME_ORIGINAL[0]
        
        algo = {"BESTA":abs(QUEUE_ABANDON_RATE_BESTA[-1]-original_abandon)}
        mse = {"BESTA":mean_squared_error(QUEUE_ABANDON_RATE_BESTA,QUEUE_ABANDON_RATE_ORIGINAL*len(QUEUE_ABANDON_RATE_BESTA))}
        rmse = {"BESTA":sqrt(mean_squared_error(QUEUE_ABANDON_RATE_BESTA,QUEUE_ABANDON_RATE_ORIGINAL*len(QUEUE_ABANDON_RATE_BESTA)))}

        for k in methods:
            if k == "perron_frobenius":
                algo["Perron Frobenius"]=abs(QUEUE_ABANDON_RATE_M1[-1]-original_abandon)
                mse["Perron Frobenius"]= mean_squared_error(QUEUE_ABANDON_RATE_M1,QUEUE_ABANDON_RATE_ORIGINAL*len(QUEUE_ABANDON_RATE_M1))
                rmse["Perron Frobenius"]= sqrt(mean_squared_error(QUEUE_ABANDON_RATE_M1,QUEUE_ABANDON_RATE_ORIGINAL*len(QUEUE_ABANDON_RATE_M1)))
            elif k == "gerschgorin":
                algo["Gerschgorin"]=abs(QUEUE_ABANDON_RATE_M2[-1]-original_abandon)
                mse["Gerschgorin"]= mean_squared_error(QUEUE_ABANDON_RATE_M2,QUEUE_ABANDON_RATE_ORIGINAL*len(QUEUE_ABANDON_RATE_M2))
                rmse["Gerschgorin"]= sqrt(mean_squared_error(QUEUE_ABANDON_RATE_M2,QUEUE_ABANDON_RATE_ORIGINAL*len(QUEUE_ABANDON_RATE_M2)))
            elif k == "partial_sum":
                algo["Partial Sum"]=abs(QUEUE_ABANDON_RATE_M3[-1]-original_abandon)
                mse["Partial Sum"]= mean_squared_error(QUEUE_ABANDON_RATE_M3,QUEUE_ABANDON_RATE_ORIGINAL*len(QUEUE_ABANDON_RATE_M3))
                rmse["Partial Sum"]= sqrt(mean_squared_error(QUEUE_ABANDON_RATE_M3,QUEUE_ABANDON_RATE_ORIGINAL*len(QUEUE_ABANDON_RATE_M3)))
            elif k == "page_rank":
                algo["Page Rank"]=abs(QUEUE_ABANDON_RATE_M4[-1]-original_abandon)
                mse["Page Rank"]= mean_squared_error(QUEUE_ABANDON_RATE_M4,QUEUE_ABANDON_RATE_ORIGINAL*len(QUEUE_ABANDON_RATE_M4))
                rmse["Page Rank"]= sqrt(mean_squared_error(QUEUE_ABANDON_RATE_M4,QUEUE_ABANDON_RATE_ORIGINAL*len(QUEUE_ABANDON_RATE_M4)))
            elif k == "chain_indexing":
                algo["Chain Indexing"]=abs(QUEUE_ABANDON_RATE_M5[-1]-original_abandon)
                mse["Chain Indexing"]= mean_squared_error(QUEUE_ABANDON_RATE_M5,QUEUE_ABANDON_RATE_ORIGINAL*len(QUEUE_ABANDON_RATE_M5))
                rmse["Chain Indexing"]= sqrt(mean_squared_error(QUEUE_ABANDON_RATE_M5,QUEUE_ABANDON_RATE_ORIGINAL*len(QUEUE_ABANDON_RATE_M5)))
            else:
                pass
        
        sorted_dict1 = OrderedDict(sorted(algo.items(), key=lambda x: x[1]))
        # sorted_dict2 = OrderedDict(sorted(mse.items(), key=lambda x: x[1]))
        # sorted_dict3 = OrderedDict(sorted(rmse.items(), key=lambda x: x[1]))
        print("Abandon Rate: ",sorted_dict1,"\n")
        # print("Abandon Rate MSE: ",sorted_dict2,"\n")
        # print("Abandon Rate RMSE: ",sorted_dict3,"\n")
        
        algo = {"BESTA":abs(QUEUE_SERVICE_RATE_BESTA[-1]-original_sr)}
        mse = {"BESTA":mean_squared_error(QUEUE_SERVICE_RATE_BESTA,QUEUE_SERVICE_RATE_ORIGINAL*len(QUEUE_SERVICE_RATE_BESTA))}
        rmse = {"BESTA":sqrt(mean_squared_error(QUEUE_SERVICE_RATE_BESTA,QUEUE_ABANDON_RATE_ORIGINAL*len(QUEUE_SERVICE_RATE_BESTA)))}
        for k in methods:
            if k == "perron_frobenius":
                algo["Perron Frobenius"]=abs(QUEUE_SERVICE_RATE_M1[-1]-original_sr)
                mse["Perron Frobenius"]= mean_squared_error(QUEUE_SERVICE_RATE_M1,QUEUE_SERVICE_RATE_ORIGINAL*len(QUEUE_SERVICE_RATE_M1))
                rmse["Perron Frobenius"]= sqrt(mean_squared_error(QUEUE_SERVICE_RATE_M1,QUEUE_SERVICE_RATE_ORIGINAL*len(QUEUE_SERVICE_RATE_M1)))
            elif k == "gerschgorin":
                algo["Gerschgorin"]=abs(QUEUE_SERVICE_RATE_M2[-1]-original_sr)
                mse["Gerschgorin"]= mean_squared_error(QUEUE_SERVICE_RATE_M2,QUEUE_SERVICE_RATE_ORIGINAL*len(QUEUE_SERVICE_RATE_M2))
                rmse["Gerschgorin"]= sqrt(mean_squared_error(QUEUE_SERVICE_RATE_M2,QUEUE_SERVICE_RATE_ORIGINAL*len(QUEUE_SERVICE_RATE_M2)))
            elif k == "partial_sum":
                algo["Partial Sum"]=abs(QUEUE_SERVICE_RATE_M3[-1]-original_sr)
                mse["Partial Sum"]= mean_squared_error(QUEUE_SERVICE_RATE_M3,QUEUE_SERVICE_RATE_ORIGINAL*len(QUEUE_SERVICE_RATE_M3))
                rmse["Partial Sum"]= sqrt(mean_squared_error(QUEUE_SERVICE_RATE_M3,QUEUE_SERVICE_RATE_ORIGINAL*len(QUEUE_SERVICE_RATE_M3)))
            elif k == "page_rank":
                algo["Page Rank"]=abs(QUEUE_SERVICE_RATE_M4[-1]-original_sr)
                mse["Page Rank"]= mean_squared_error(QUEUE_SERVICE_RATE_M4,QUEUE_SERVICE_RATE_ORIGINAL*len(QUEUE_SERVICE_RATE_M4))
                rmse["Page Rank"]= sqrt(mean_squared_error(QUEUE_SERVICE_RATE_M4,QUEUE_SERVICE_RATE_ORIGINAL*len(QUEUE_SERVICE_RATE_M4)))
            elif k == "chain_indexing":
                algo["Chain Indexing"]=abs(QUEUE_SERVICE_RATE_M5[-1]-original_sr)
                mse["Chain Indexing"]= mean_squared_error(QUEUE_SERVICE_RATE_M5,QUEUE_SERVICE_RATE_ORIGINAL*len(QUEUE_SERVICE_RATE_M5))
                rmse["Chain Indexing"]= sqrt(mean_squared_error(QUEUE_SERVICE_RATE_M5,QUEUE_SERVICE_RATE_ORIGINAL*len(QUEUE_SERVICE_RATE_M5)))
            else:
                pass
                
        sorted_dict1 = OrderedDict(sorted(algo.items(), key=lambda x: x[1]))
        # sorted_dict2 = OrderedDict(sorted(mse.items(), key=lambda x: x[1]))
        # sorted_dict3 = OrderedDict(sorted(rmse.items(), key=lambda x: x[1]))
        print("Service Rate: ",sorted_dict1,"\n")
        # print("Service Rate MSE: ",sorted_dict2,"\n")
        # print("Service Rate RMSE: ",sorted_dict3,"\n")

        algo = {"BESTA":abs(QUEUE_SERVICE_TIME_BESTA[-1]-original_st)}
        mse["BESTA"]= mean_squared_error(QUEUE_SERVICE_TIME_BESTA,QUEUE_SERVICE_TIME_ORIGINAL*len(QUEUE_SERVICE_TIME_BESTA))
        rmse["BESTA"]= sqrt(mean_squared_error(QUEUE_SERVICE_TIME_M5,QUEUE_SERVICE_TIME_ORIGINAL*len(QUEUE_SERVICE_TIME_BESTA)))
        
        for k in methods:
            if k == "perron_frobenius":
                algo["Perron Frobenius"]=abs(QUEUE_SERVICE_TIME_M1[-1]-original_st)
                mse["Perron Frobenius"]= mean_squared_error(QUEUE_SERVICE_TIME_M1,QUEUE_SERVICE_TIME_ORIGINAL*len(QUEUE_SERVICE_TIME_M1))
                rmse["Perron Frobenius"]= sqrt(mean_squared_error(QUEUE_SERVICE_TIME_M1,QUEUE_SERVICE_TIME_ORIGINAL*len(QUEUE_SERVICE_TIME_M1)))
            elif k == "gerschgorin":
                algo["Gerschgorin"]=abs(QUEUE_SERVICE_TIME_M2[-1]-original_st)
                mse["Gerschgorin"]= mean_squared_error(QUEUE_SERVICE_TIME_M2,QUEUE_SERVICE_TIME_ORIGINAL*len(QUEUE_SERVICE_TIME_M2))
                rmse["Gerschgorin"]= sqrt(mean_squared_error(QUEUE_SERVICE_TIME_M2,QUEUE_SERVICE_TIME_ORIGINAL*len(QUEUE_SERVICE_TIME_M2)))
            elif k == "partial_sum":
                algo["Partial Sum"]=abs(QUEUE_SERVICE_TIME_M3[-1]-original_st)
                mse["Partial Sum"]= mean_squared_error(QUEUE_SERVICE_TIME_M3,QUEUE_SERVICE_TIME_ORIGINAL*len(QUEUE_SERVICE_TIME_M3))
                rmse["Partial Sum"]= sqrt(mean_squared_error(QUEUE_SERVICE_TIME_M3,QUEUE_SERVICE_TIME_ORIGINAL*len(QUEUE_SERVICE_TIME_M3)))
            elif k == "page_rank":
                algo["Page Rank"]=abs(QUEUE_SERVICE_TIME_M4[-1]-original_st)
                mse["Page Rank"]= mean_squared_error(QUEUE_SERVICE_TIME_M4,QUEUE_SERVICE_TIME_ORIGINAL*len(QUEUE_SERVICE_TIME_M4))
                rmse["Page Rank"]= sqrt(mean_squared_error(QUEUE_SERVICE_TIME_M4,QUEUE_SERVICE_TIME_ORIGINAL*len(QUEUE_SERVICE_TIME_M4)))
            elif k == "chain_indexing":
                algo["Chain Indexing"]=abs(QUEUE_SERVICE_TIME_M5[-1]-original_st)
                mse["Chain Indexing"]= mean_squared_error(QUEUE_SERVICE_TIME_M5,QUEUE_SERVICE_TIME_ORIGINAL*len(QUEUE_SERVICE_TIME_M5))
                rmse["Chain Indexing"]= sqrt(mean_squared_error(QUEUE_SERVICE_TIME_M5,QUEUE_SERVICE_TIME_ORIGINAL*len(QUEUE_SERVICE_TIME_M5)))
            else:
                pass

        sorted_dict1 = OrderedDict(sorted(algo.items(), key=lambda x: x[1]))
        # sorted_dict2 = OrderedDict(sorted(mse.items(), key=lambda x: x[1]))
        # sorted_dict3 = OrderedDict(sorted(rmse.items(), key=lambda x: x[1]))
        print("Service Time: ",sorted_dict1,"\n")
        # print("Service Time MSE: ",sorted_dict2,"\n")
        # print("Service Time RMSE: ",sorted_dict3,"\n")
        
        algo = {"BESTA":abs(QUEUE_WAITING_TIME_BESTA[-1]-original_wt)}
        mse["BESTA"]= mean_squared_error(QUEUE_WAITING_TIME_BESTA,QUEUE_WAITING_TIME_ORIGINAL*len(QUEUE_WAITING_TIME_BESTA))
        rmse["BESTA"]= sqrt(mean_squared_error(QUEUE_WAITING_TIME_BESTA,QUEUE_WAITING_TIME_ORIGINAL*len(QUEUE_WAITING_TIME_BESTA)))
        for k in methods:
            if k == "perron_frobenius":
                algo["Perron Frobenius"]=abs(QUEUE_WAITING_TIME_M1[-1]-original_wt)
                mse["Perron Frobenius"]= mean_squared_error(QUEUE_WAITING_TIME_M1,QUEUE_WAITING_TIME_ORIGINAL*len(QUEUE_WAITING_TIME_M1))
                rmse["Perron Frobenius"]= sqrt(mean_squared_error(QUEUE_WAITING_TIME_M1,QUEUE_WAITING_TIME_ORIGINAL*len(QUEUE_WAITING_TIME_M1)))
            elif k == "gerschgorin":
                algo["Gerschgorin"]=abs(QUEUE_WAITING_TIME_M2[-1]-original_wt)
                mse["Gerschgorin"]= mean_squared_error(QUEUE_WAITING_TIME_M2,QUEUE_WAITING_TIME_ORIGINAL*len(QUEUE_WAITING_TIME_M2))
                rmse["Gerschgorin"]= sqrt(mean_squared_error(QUEUE_WAITING_TIME_M2,QUEUE_WAITING_TIME_ORIGINAL*len(QUEUE_WAITING_TIME_M2)))
            elif k == "partial_sum":
                algo["Partial Sum"]=abs(QUEUE_WAITING_TIME_M3[-1]-original_wt)
                mse["Partial Sum"]= mean_squared_error(QUEUE_WAITING_TIME_M3,QUEUE_WAITING_TIME_ORIGINAL*len(QUEUE_WAITING_TIME_M3))
                rmse["Partial Sum"]= sqrt(mean_squared_error(QUEUE_WAITING_TIME_M3,QUEUE_WAITING_TIME_ORIGINAL*len(QUEUE_WAITING_TIME_M3)))
            elif k == "page_rank":
                algo["Page Rank"]=abs(QUEUE_WAITING_TIME_M4[-1]-original_wt)
                mse["Page Rank"]= mean_squared_error(QUEUE_WAITING_TIME_M4,QUEUE_WAITING_TIME_ORIGINAL*len(QUEUE_WAITING_TIME_M4))
                rmse["Page Rank"]= sqrt(mean_squared_error(QUEUE_WAITING_TIME_M4,QUEUE_WAITING_TIME_ORIGINAL*len(QUEUE_WAITING_TIME_M4)))
            elif k == "chain_indexing":
                algo["Chain Indexing"]=abs(QUEUE_WAITING_TIME_M5[-1]-original_wt)
                mse["Chain Indexing"]= mean_squared_error(QUEUE_WAITING_TIME_M5,QUEUE_WAITING_TIME_ORIGINAL*len(QUEUE_WAITING_TIME_M5))
                rmse["Chain Indexing"]= sqrt(mean_squared_error(QUEUE_WAITING_TIME_M5,QUEUE_WAITING_TIME_ORIGINAL*len(QUEUE_WAITING_TIME_M5)))
            else:
                pass

        sorted_dict1 = OrderedDict(sorted(algo.items(), key=lambda x: x[1]))
        # sorted_dict2 = OrderedDict(sorted(mse.items(), key=lambda x: x[1]))
        # sorted_dict3 = OrderedDict(sorted(rmse.items(), key=lambda x: x[1]))
        print("Waiting Time: ",sorted_dict1,"\n")
        # print("Waiting Time MSE: ",sorted_dict2,"\n")
        # print("Waiting Time RMSE: ",sorted_dict3,"\n")