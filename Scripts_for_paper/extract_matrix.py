### to generate the ergodic matrix into the file day1_day2.dat from the file sended by Yuxin Scenario2_Last_5_Days.dat
### >python extract_matrix.py Scenario2_Last_5_Days.dat 1 2
### to generate the ergodic matrix into the file day2_day3.dat from the file sended by Yuxin Scenario2_Last_5_Days.dat
### >python extract_matrix.py Scenario2_Last_5_Days.dat 2 3
### etc...

import numpy as np
import sys, os
import seaborn as sns
import matplotlib.pyplot as plt
from pprint import pprint
from itertools import product
from sinkhorn_knopp import sinkhorn_knopp as skp

def generate_combinations(n, m):
    """Génère toutes les combinaisons possibles de n éléments prenant des valeurs de 0 à m."""
    values = range(m + 1)
    return [a for a in product(values, repeat=n) if sum(a) <= m]

def initialize_matrix_with_small_values(n, scale=1e-6):
    """Initialise une matrice carrée de taille (m+1)^n avec des petites valeurs proches de 0."""
    # Utilise une distribution normale centrée sur 0 avec une petite variance
    return np.random.uniform(low=0.0, high=scale, size=(n,n))

def check_row_sums(matrix, tolerance=1e-6):
    """
    Vérifie que la somme des lignes de la matrice est égale à 1.

    Args:
        matrix (numpy.ndarray): La matrice à vérifier.
        tolerance (float): La tolérance pour la vérification de l'égalité (par défaut 1e-6).

    Returns:
        bool: True si la somme des lignes est égale à 1 dans la tolérance spécifiée, sinon False.
    """
    row_sums = np.sum(matrix, axis=1)
    return np.all(np.abs(row_sums - 1) < tolerance)

def check_column_sums(matrix, tolerance=1e-2):
    """
    Vérifie que la somme des colonnes de la matrice est égale à 1.

    Args:
        matrix (numpy.ndarray): La matrice à vérifier.
        tolerance (float): La tolérance pour la vérification de l'égalité (par défaut 1e-6).

    Returns:
        bool: True si la somme des colonnes est égale à 1 dans la tolérance spécifiée, sinon False.
    """
    column_sums = np.sum(matrix, axis=0)
    return np.all(np.abs(column_sums - 1) < tolerance)

def read_yuxin_file(filename,day1,day2):
    L = []

    with open(filename,'r') as f:
        for line in f.readlines():
            
            l,prob = line.split(') ')
    
            cash = float(l.split(',')[-1])
            prob = prob.replace("\n","")

            l = eval(l+')')
            if l[0]==day1:
                L.append([f"({l[1]},{l[2]},{l[3]})"])
            elif l[0]==day2:
                L[-1].append(f"({l[1]},{l[2]},{l[3]});{prob}")
            
    D = {}

    ## update of value given by Yuxin
    for sub_list in L:
        start_stat = sub_list[0]
        for info in sub_list[1:]:
            end_state, prob = info.split(';')
            if (start_stat, end_state) not in D.keys(): 
                D[(start_stat, end_state)] = eval(prob)
            else:
                ### add prob
                D[(start_stat, end_state)] += eval(prob)

            # print(f"{start_stat} {end_state} {prob}")

    return D

def write_matrix_to_file(matrix, combinations, filename):
    """
    """
    
    with open(filename, 'w') as f:
        for i,s1 in enumerate(combinations):
            for j,s2 in enumerate(combinations):
                value = matrix[i, j]
                # Formater les indices et la valeur
                line = f"{''.join(map(str, s1))} {''.join(map(str, s2))} {value:.10f}\n"
                f.write(line)

def afficher_carte_de_chaleur(matrice, row_labels, col_labels):
    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(matrice, annot=False, fmt=".2e", cmap='viridis', xticklabels = row_labels, yticklabels= col_labels)
    ax.xaxis.tick_top()
    # plt.title("Ergotic Matrix")
    plt.xlabel("State")
    plt.ylabel("State")
    plt.show()

def print_help():
    help_text = """
    Usage: python extract_matrix.py [options] <source_filename> <day1> <day2>

    Options:
    -H, --help        : Displays this help message and exits the program.

    Arguments:
    <source_filename> : The source file to process. This file is provided by the Yuxoin program.
    <day1>            : The day n for processing, an integer.
    <day2>            : The day n+1 for processing, an integer.

    Description:
    This script takes three arguments:
    1. The source file name.
    2. Day n.
    3. Day n+1.

    The output file name will be generated based on the specified days, with the following format:
    day{day1}_day{day2}.dat

    Example:
    python script.py input_file.dat 1 2
    This will generate an output file named day1_day2.dat.
    """
    print(help_text)

def main():

    ### How to use
    ### generate matrix day3_day4.dat for day 3 to day 4
    ### >python extract_matrix.py Scenario2_Last_5_Days.dat 3 4

    if len(sys.argv) < 2:
        print("Error: Incorrect number of arguments.")
        print_help()
        sys.exit(1)

    if sys.argv[1] in ['-h', '--help']:
        print_help()
        sys.exit(0)

    if len(sys.argv) != 4:
        print("Error: Incorrect number of arguments.")
        print_help()
        sys.exit(1)

    try:
        source_filename = sys.argv[1]  # Source file from Yuxoin prog
        
        if not os.path.exists(source_filename):
            print("Error: First param must be a file")
            print_help()
            sys.exit(1)

        day1 = int(sys.argv[2])       # Day n
        day2 = int(sys.argv[3])       # Day n+1

        # Construct the full path to the output file
        current_dir = os.path.dirname(os.path.abspath(__file__))
        yuxin_dir = os.path.join(current_dir, '..', 'Matrix', 'Yuxin')
        yuxin_matrix_dir = yuxin_dir if os.path.exists(yuxin_dir) else ""
        
        output_fn = os.path.join(yuxin_matrix_dir,f"{os.path.splitext(source_filename)[0]}_day{day1}_to_day{day2}.dat")  # Output file name

    except ValueError:
        print("Error: Arguments <day1> and <day2> must be integers.")
        print_help()
        sys.exit(1)
    
    n = 3  # Number of stocks
    m = day2  # max val for each stock
    
    # Générer toutes les combinaisons possibles
    combinations = generate_combinations(n, m)
    
    matrix_size = len(combinations)

    # Initialiser la matrice avec des petites valeurs proches de 0
    matrix = initialize_matrix_with_small_values(matrix_size)
    
    print(f"Size of the matrix for day {day1} to day {day2} : {matrix_size}x{matrix_size}")

    # Affichage des combinaisons pour vérification
    # print("Toutes les combinaisons possibles :")
    # for idx, comb in enumerate(combinations):
        # print(f"Index {idx} : {comb}")
    
    # Exemple de manipulation : définir une valeur de transition
    # index_from = combinations.index((2, 1, 0))
    # index_to = combinations.index((0, 0, 1))
    # matrix[index_from, index_to] = 1.0

    print(f"Extraction of file {source_filename}...")
    D = read_yuxin_file(source_filename, day1, day2)
    pprint(sorted(D.items(), key=lambda item:item[1], reverse=True))

    print("Done!")

    ### update the matrix from D
    for k,v in D.items():
        index_from = combinations.index(eval(k[0]))
        index_to = combinations.index(eval(k[1]))
        matrix[index_from, index_to] = v
   
    ### make matrix ergodic
    sk = skp.SinkhornKnopp()
    P_ds = sk.fit(matrix)

    # Régler les options d'affichage pour une meilleure lisibilité
    # np.set_printoptions(precision=4, suppress=True)

    # Afficher la matrice transformée
    # print(P_ds)

    # Vérifier que la somme des lignes et colonne fait 1
    assert check_row_sums(P_ds)  
    assert check_column_sums(P_ds)
    
    ### print map of matrix
    labels = ["".join(map(str,a)) for a in combinations]
    afficher_carte_de_chaleur(P_ds, labels, labels)

    ### write output file
    write_matrix_to_file(P_ds, combinations, output_fn)
    
    print("List of considered states:")
    print(combinations)

if __name__ == "__main__":
    main()