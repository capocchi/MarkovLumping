import numpy as np

class Queue:
    
    def __init__(self, arrival_rate, service_rate, capacity):
        """
        arrival_rate : le taux d'arrivée des clients (en clients/seconde)
        service_rate : le taux de service des clients (en clients/seconde)
        capacity : la capacité maximale de la file d'attente
        """
        self.arrival_rate = arrival_rate
        self.service_rate = service_rate
        self.capacity = capacity
        self.state = 0

    def transition_matrix(self):
        """
        génère la matrice de transition ergodique de la file d'attente. 
        Cette matrice modélise le comportement de la file d'attente et permet de calculer des statistiques telles que le temps moyen d'attente.
        """
        transition = np.zeros((self.capacity + 1, self.capacity + 1))
        for i in range(self.capacity):
            transition[i, i+1] = self.arrival_rate
            transition[i+1, i] = min(i, self.service_rate)
            transition[i, i] = self.arrival_rate + max(0, self.service_rate - i)
        transition[self.capacity, self.capacity] = self.arrival_rate + self.service_rate
        return transition

if __name__ == '__main__':
    """
    Un exemple d'utilisation de cette classe serait de modéliser une file 
    d'attente dans un centre d'appels avec une capacité maximale de 10 appels, 
    un taux d'arrivée de 5 appels par minute et un taux de service de 3 appels par minute:
    """
    queue = Queue(arrival_rate=5/60, service_rate=3/60, capacity=10)
    transition_matrix = queue.transition_matrix()
    print(transition_matrix)
