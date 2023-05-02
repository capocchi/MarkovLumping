#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This script is used to compute the best partition of a (pykov) ergotic Markov chain using the KL rate. See the __main__ for the usage.
# Copyright (C) 2021  Laurent Capocchi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Email: capocchi@univ-corse.fr

#######################################################################
### Name: queuing.py
### Author: L. Capocchi
### Version: 1.0
### Description: script to reduce a queue depending ion the the main_loop algo. See the __main__ for the usage. 
### Dependencies: pykov, networkx, numpy, matplotlib
### Python version: 3.9
### Date: 28/02/2023
#######################################################################


import numpy as np
import networkx as nx
from sklearn import preprocessing

from main_loop import *

PLOT = True

def generate_queue_markov_chain(arrival_rate, service_rate, num_servers):
    """
    Génère une chaîne de Markov correspondant à une file d'attente
    avec un taux d'arrivée, un taux de service et un nombre de serveurs donnés.
    """
    # Calcul du nombre total d'états
    num_states = num_servers + 1
    
    # Initialisation de la matrice de transition
    transition_matrix = np.zeros((num_states, num_states))
    
    # Calcul des taux de transition
    for i in range(num_states):
        for j in range(num_states):
            if i < j:
                transition_matrix[i, j] = arrival_rate * (num_servers - i)
            elif i == j:
                transition_matrix[i, j] = arrival_rate * (num_servers - i) + service_rate * i
            else:
                transition_matrix[i, j] = service_rate * i
                
    # Normalisation de la matrice de transition
    row_sums = np.sum(transition_matrix, axis=1)
    transition_matrix = np.divide(transition_matrix, row_sums[:, np.newaxis], where=row_sums[:, np.newaxis]!=0)
    
    return transition_matrix

def generate_mm1_queue_chain(lambda_val, mu_val, num_states):
    """
    Génère une chaîne de Markov pour une file d'attente M/M/1 avec un taux d'arrivée lambda_val, un taux de service mu_val,
    et num_states états. Chaque état représente le nombre de clients dans la file d'attente et dans le système de service.
    Cette fonction prend trois arguments en entrée : le taux d'arrivée des clients lambda_val, le taux de service des clients mu_val,
     et le nombre d'états dans la chaîne de Markov num_states. La fonction renvoie une matrice de transition de taille (num_states, num_states) 
     qui représente la chaîne de Markov pour une file d'attente M/M/1. Chaque état de la chaîne de Markov représente le nombre de clients dans la 
     file d'attente et dans le système de service.

    Args:
    - lambda_val (float): Taux d'arrivée des clients
    - mu_val (float): Taux de service des clients
    - num_states (int): Nombre d'états dans la chaîne de Markov
    
    Returns:
    - np.array: Une matrice de transition de taille (num_states, num_states) qui représente la chaîne de Markov.
    """
    arrival_rates = np.zeros(num_states) # taux d'arrivée dans chaque état
    service_rates = np.zeros(num_states) # taux de service dans chaque état

    # Calcul du taux d'arrivée et du taux de service pour chaque état
    for i in range(num_states):
        arrival_rates[i] = lambda_val
        service_rates[i] = min(i, 1) * mu_val
    
    # Initialisation de la matrice de transition
    transition_matrix = np.zeros((num_states, num_states))

    # Calcul des probabilités de transition entre les états
    for i in range(num_states):
        transition_matrix[i, i] = sum(arrival_rates[:i]) + sum(service_rates[i+1:])
        if i < num_states - 1:
            transition_matrix[i, i+1] = arrival_rates[i]
        if i > 0:
            transition_matrix[i, i-1] = service_rates[i]

    # Normalisation des probabilités de transition
    for i in range(num_states):
        transition_matrix[i, :] /= sum(transition_matrix[i, :])

    return transition_matrix

def calculate_queue_utilization(transition_matrix):
    """
    Calcule le facteur de charge d'une file d'attente
    à partir de sa matrice de transition.
    """

    num_servers = transition_matrix.shape[0] - 1
    service_rate = 0
    
    # Calcul du taux de service moyen
    for i in range(num_servers + 1):
        for j in range(num_servers + 1):
            if i > 0 and i == j:
                service_rate += transition_matrix[i, j]
    
    # Calcul du facteur de charge
    arrival_rate = transition_matrix[0, 1:].sum()
    utilization_factor = arrival_rate / (num_servers * service_rate)
    
    return utilization_factor

def waiting_time_from_transition_matrix(transition_matrix):
    """
    Calcule le temps moyen d'attente dans une file d'attente à partir de sa matrice de transition.

    Args:
        Q (np.ndarray): La matrice de transition de la file d'attente.

    Returns:
        float: Le temps moyen d'attente dans la file d'attente.
    """
    
    # On vérifie si la matrice de transition est bien stochastique
    #if not np.all(np.linalg.eigvals(transition_matrix) > 0):
    #    raise ValueError("La matrice de transition n'est pas stochastique.")
    
    # On calcule la matrice Q, la matrice de transition sans la diagonale
    Q = transition_matrix[:-1,:-1]
    # On calcule la matrice R, la dernière colonne de la matrice de transition
    R = transition_matrix[:-1,-1]
    # On calcule la matrice I, la matrice identité de la même taille que Q
    I = np.identity(Q.shape[0])
    # On calcule la matrice N, la matrice fondamentale
    N = np.linalg.inv(I - Q)
    # On calcule le vecteur de temps moyen d'attente
    T = np.dot(N, R)
    # On retourne le temps moyen d'attente pour chaque état
    return T

def transform(e):
    if '/' in e:
        e = e.split('/')[0]

    return int(e)

def chainToNPArray(P):
    """ Transfrom the pykov chain into np array
    """
    
    if type(P) is dict :
        PP = {}
        for k,v in P.items():
            PP[tuple(map(lambda a: a.replace('NS',''),k))] = v 
        P = pykov.Chain(PP)

    D = {}
    for i,partition in enumerate(sorted(list(P.states()), key=transform)):
        D[partition] = i

    # final matrix to return P.states()xP.states()
    L = np.array([[0]*len(P.states())]*len(P.states()),dtype=float)

    # make the final matrix L from D
    for k in list(P.keys()):
        a,b = k
        L[int(D[a]),int(D[b])] = P[a,b]

    return L



def getStat(P, p=False):
    """
    Dans ce code, nous avons défini la matrice de transition égotique P de la file d'attente. 
    Nous avons ensuite utilisé la formule de Little pour calculer le temps d'attente moyen en utilisant 
    la matrice d'occupation M. Le temps de service moyen a été calculé en utilisant la formule de la distribution 
    de probabilité de la variable aléatoire exponentielle. Le taux de service a été calculé en utilisant la probabilité 
    de transition de l'état 1 vers l'état 2, qui représente l'état dans lequel il n'y a pas d'éléments en attente.
     Le taux d'abandon a été calculé en utilisant la probabilité de transition de l'état 1 vers l'état 4, 
     qui représente l'état dans lequel les éléments quittent la file d'attente sans être servis.
      Les résultats sont affichés à l'écran. Notez que ces calculs dépendent de la structure de la matrice de transition égotique, 
      et doivent être ajustés en conséquence pour d'autres types de file d'attente.
    """
    try:
        # Calcul du temps d'attente moyen
        I = np.identity(len(P))
        M = np.linalg.inv(I - P)
        avg_waiting_time = M.sum(axis=1)[1] / M[1, 1]
    except:
        avg_waiting_time = -1

    try:
        # Calcul du temps de service moyen
        avg_service_time = 1 / (1 - P[0, 0])
    except:
        avg_service_time = 10e10

    # Calcul du taux de service
    service_rate = 1 - P[0, 0]

    # Calcul du taux d'abandon
    abandon_rate = P[0, -1]

    # Affichage des résultats
    if p:
        print("Temps d'attente moyen : {:.4f}".format(avg_waiting_time))
        print("Temps de service moyen : {:.4f}".format(avg_service_time))
        print("Taux de service : {:.4f}".format(service_rate))
        print("Taux d'abandon : {:.4f}".format(abandon_rate))

    return (avg_waiting_time,avg_service_time,service_rate,abandon_rate) 

def getStat2(transition_matrix):
    # extraire le nombre d'états du système
    num_states = transition_matrix.shape[0]

    # calculer la probabilité stationnaire du système en résolvant
    # l'équation pi * P = pi, où P est la matrice de transition
    # et pi est la distribution de probabilité stationnaire
    pi = np.linalg.solve(np.transpose(transition_matrix) - np.eye(num_states), np.zeros(num_states) + 1)
    
    # extraire les taux d'arrivée et de service à partir de la matrice de transition
    arrival_rate = transition_matrix[0, 1]
    service_rate = -np.diag(transition_matrix)[0]
    
    # calculer le temps d'attente moyen et le temps de service moyen
    average_waiting_time = pi[1] / (arrival_rate * (1 - pi[0]))
    average_service_time = 1 / service_rate
    
    # calculer le taux de service et le taux d'abandon
    service_rate_percent = service_rate / arrival_rate * 100
    abandonment_rate_percent = pi[0] * 100

        # Affichage des résultats
    print("Temps d'attente moyen : {:.2f}".format(average_waiting_time))
    print("Temps de service moyen : {:.2f}".format(average_service_time))
    print("Taux de service : {:.2f}".format(service_rate_percent))
    print("Taux d'abandon : {:.2f}".format(abandonment_rate_percent))
    
    # retourner les résultats
    return (average_waiting_time, average_service_time, service_rate_percent, abandonment_rate_percent)

if __name__ == '__main__':

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
        # starting time
        start1 = time.time()
        
        num_servers = int(n)
        service_rate = float(sr)
        arrival_rate = float(ar)

        # Le temps moyen d'attente dans la file d'attente
        #mean_waiting_time = (1 / (service_rate - arrival_rate)) * (1 - (num_servers * (arrival_rate / service_rate))**num_servers * (1 - (arrival_rate / service_rate))**(num_servers + 1))

        # facteur de charge 
        #rho = arrival_rate / (num_servers * service_rate)

        queue_transition_matrix = generate_mm1_queue_chain(arrival_rate,service_rate,num_servers)

        #queue_transition_matrix= generate_queue_markov_chain(arrival_rate,service_rate,num_servers)
        #queue_transition_matrix = Queue(arrival_rate, service_rate, num_servers).transition_matrix()
        
        P = pykov.Chain()
        for i,row in enumerate(queue_transition_matrix):
            for j,val in enumerate(row):
                P[(str(i),str(j))]=val
        
        S = tuple(sorted(list(P.states())))
        kl,p,Q = next(getMFTPAnalysis2(S,P))
        n = len(S)

        ### transform Chain to matrix
        A = chainToNPArray(P)

        ### get queue stat
        avg_waiting_time,avg_service_time,service_rate,abandon_rate = getStat(A)

        #A_reduced_by_perron_frobenius = perron_frobenius_reduction(A)
        #P_reduced_by_perron_frobenius = pykov.Chain()
        #for i,row in enumerate(A_reduced_by_perron_frobenius):
        #    for j,val in enumerate(row):
        #        P_reduced_by_perron_frobenius[(str(i),str(j))]=val
 
        #S_reduced_by_perron_frobenius = tuple(sorted(list(P_reduced_by_perron_frobenius.states())))
        #kl_P_reduced_by_perron_frobenius,_,_ = next(getMFTPAnalysis2(S_reduced_by_perron_frobenius,P_reduced_by_perron_frobenius))
        #print(kl_P_reduced_by_perron_frobenius)

        #queue_utilization = calculate_queue_utilization(A)
        #queue_waiting_time = waiting_time_from_transition_matrix(A)

        #print(f"initial queue utilization from P/formula: {queue_utilization}/{rho}")
        #print(f"initial mean waiting time from P/formula: {queue_waiting_time}/{mean_waiting_time}")

        # starting time
        end1 = time.time()

        # total time taken
        print(f"\nRuntime of the first phase of the algorithm BESTA is {end1 - start1}s")

        if PLOT:
            X = []
            K_L = []
            QUEUE_WAITING_TIME = []
            QUEUE_SERVICE_TIME = []
            QUEUE_SERVICE_RATE = []
            QUEUE_ABANDON_RATE = []
        
        if STAT:    
            STEADY = {n:P.steady()}
            if not PLOT:
                K_L = []
                
        if PLOT:
            X.append(n)
            K_L.append(kl)
            QUEUE_WAITING_TIME.append(avg_waiting_time)
            QUEUE_SERVICE_TIME.append(avg_service_time)
            QUEUE_SERVICE_RATE.append(service_rate)
            QUEUE_ABANDON_RATE.append(abandon_rate)

        ### condition for table 2
        cond = "n>2"

        ### condition for table3
        #kl=new_kl=diff=0.0
        #cond = "new_kl <= kl*(1+0.5) or kl==0.0"

        #print(np.array_equal(L, generate_queue_markov_chain(arrival_rate,service_rate,num_servers)))
        
        ### stopping condition
        while(eval(cond)):
            
            # starting time
            start2 = time.time()

            ### P must be ergotic i.e. the transition matrix must be irreducible and acyclic.            
            G = nx.DiGraph(list(P.keys()), directed=True)
            nx.strongly_connected_components(G)
            assert nx.is_strongly_connected(G) and nx.is_aperiodic(G), f"Matrix is not ergotic!"
            
            ### just for stdout
            d = trace(n,kl,p)
            
            ### P transition matrix of the Markoc chain to lump
            P = pykov.Chain()
            for k,v in Q.items():
                P[(d[k[0]],d[k[1]])] = v

            A = chainToNPArray(P)
        
            try:
                avg_waiting_time,avg_service_time,service_rate,abandon_rate = getStat(A)

            except Exception as e:
                print("Une exception s'est produite :", e)

            ### set of states
            S = tuple(sorted(P.states()))

            if "n>2" not in cond:
                kl = new_kl
                
            ### mfpt analisys (REDUCED function that call FOP)
            new_kl,p,Q = next(getMFTPAnalysis2(S,P))
            
            ### update variables
            n = len(S)
           
            if "n>2" in cond:
                kl = new_kl

            #print(new_kl, kl*(1+0.5))
            
            # starting time
            end2 = time.time()

            # total time taken
            print(f"\nRuntime of the 'while' phase of the algorithm BESTA for n={n} is {end2 - start2}s")

            if PLOT:
                X.append(n)
                K_L.append(kl)
                QUEUE_WAITING_TIME.append(avg_waiting_time)
                QUEUE_SERVICE_TIME.append(avg_service_time)
                QUEUE_SERVICE_RATE.append(service_rate)
                QUEUE_ABANDON_RATE.append(abandon_rate)
                #displayGraph(dict(P))
            
            if STAT:
                STEADY[n]=P.steady()
                if not PLOT:
                    K_L.append(kl)
                
            if WRITE_FILE:
                fn = os.path.join(os.pardir,'Matrix',f"{n}x{n}.dat")
                if os.path.exists(fn):
                    os.remove(fn)
                f = open(fn,'w')
                for k,v in dict(P).items():
                    f.write(f"{k[0]} {k[1]} {v} \n")
                f.close()
         
        ### just for stdout
        trace(n,kl,p)
        A = chainToNPArray(P)
        #print(f"queue utilization: {calculate_queue_utilization(A)}")
        #print(f"mean waiting time: {waiting_time_from_transition_matrix(A)}")
        
        try:
            avg_waiting_time,avg_service_time,service_rate,abandon_rate = getStat(A)
        except:
            pass

        # end time
        end3 = time.time()

        # total time taken
        print(f"\nRuntime of the program is {end3 - start1}s")
     
        if STAT:    
            import pandas as pd
            import pprint
            s = pd.Series(K_L)
            print(s.describe())
            pprint.pprint(STEADY)
            
        if PLOT:
            X = list(map(lambda a : f"{a} x {a}", map(str,X)))

            normalized_kl = preprocessing.normalize([K_L])
            normalized_QUEUE_WAITING_TIME = preprocessing.normalize([QUEUE_WAITING_TIME])
            normalized_QUEUE_SERVICE_TIME = preprocessing.normalize([QUEUE_SERVICE_TIME])
            normalized_QUEUE_SERVICE_RATE = preprocessing.normalize([QUEUE_SERVICE_RATE])
            normalized_QUEUE_ABANDON_RATE = preprocessing.normalize([QUEUE_ABANDON_RATE])

            # # Créer une figure avec deux sous-graphiques
            # fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(nrows=5, ncols=1)

            # # Tracer la première courbe sur le premier sous-graphique     
            # ax1.plot(X, normalized_kl[0], label="KL",marker="o")
            # ax1.set_title('KL')
            # ax1.tick_params(rotation=45)
            # ax1.grid()

            # # Tracer la deuxième courbe sur le deuxième sous-graphique
            # ax2.plot(X, normalized_QUEUE_WAITING_TIME[0], label="Queue waiting time",marker="o")
            # ax2.set_title('Queue waiting time')
            # ax2.tick_params(rotation=45)
            # ax2.grid()

            # # Tracer la troisième courbe sur le deuxième sous-graphique
            # ax3.plot(X, normalized_QUEUE_SERVICE_TIME[0], label="Queue service time",marker="o")
            # ax3.set_title('Queue service time')
            # ax3.tick_params(rotation=45)
            # ax3.grid()

            #  # Tracer la quatrième courbe sur le deuxième sous-graphique
            # ax4.plot(X, normalized_QUEUE_SERVICE_RATE[0], label="Queue service rate",marker="o")
            # ax4.set_title('Queue service rate')
            # ax4.tick_params(rotation=45)
            # ax4.grid()

            #  # Tracer la cinquièlme courbe sur le deuxième sous-graphique
            # ax5.plot(X, normalized_QUEUE_ABANDON_RATE[0], label="Queue abandon rate",marker="o")
            # ax5.set_title('Queue abandon rate')
            # ax5.tick_params(rotation=45)
            # ax5.grid()

            # # Introduire un espace entre les subplots
            # fig.subplots_adjust(hspace=1)

            if STAT:
                plt.plot(X,[s.describe()['std']]*len(X), label='std', linestyle="-")

            plt.plot(X,normalized_kl[0], label="KL", marker="o")
            plt.plot(X,normalized_QUEUE_WAITING_TIME[0], label="Queue waiting time", marker=".")
            plt.plot(X,normalized_QUEUE_SERVICE_TIME[0], label="Queue service time", marker="1")
            plt.plot(X,normalized_QUEUE_SERVICE_RATE[0], label="Queue service rate", marker="2")
            plt.plot(X,normalized_QUEUE_ABANDON_RATE[0], label="Queue abandon rate", marker="3")

            plt.xticks(rotation=45)
            
            # show a legend on the plot
            plt.legend()
            #plt.axis([max(X),min(X),min(normalized_kl[0]),max(normalized_kl[0])])
            plt.grid()
            
            plt.ylabel("Normalized Amplitude")
            plt.xlabel("Matrix NxN")

            # function to show the plot
            plt.show()
    