# combined_script.py

# Importer les bibliothèques nécessaires
import os
import sys
import networkx as nx

# --- Extraction Matrix (Depuis extract_matrix.py) ---
def extract_matrix(input_file, start_day, end_day):
    # Logique complète de `extract_matrix.py` ici
    # ...
    print(f"Matrix extracted for days {start_day} to {end_day}")
    # Retourner la sortie nécessaire, comme le nom du fichier généré
    output_file = f"../Matrix/Yuxin/Scenario1_First_{end_day}_Days_day{start_day}_to_day{end_day}.dat"
    return output_file

# --- Main Loop (Depuis main_loop.py) ---
def main_loop(matrix_file):
    # Logique complète de `main_loop.py` ici
    # ...
    print(f"Main loop processing done for matrix file: {matrix_file}")
    # Retourner le nom du fichier exporté par exemple
    exported_file = "graph_output.graphml"  # Exemple
    return exported_file

# --- Fonction de fusion ---
def execute_combined_process(end_depth):
    # Exécuter l'extraction de la matrice
    input_file = f"Scenario1_First_{end_depth}_Days.dat"
    matrix_file = extract_matrix(input_file, end_depth - 1, end_depth)
    
    # Exécuter la boucle principale avec le fichier de matrice généré
    exported_file = main_loop(matrix_file)
    
    return exported_file

# --- Appel principal ---
if __name__ == "__main__":
    end_depth = int(sys.argv[1])  # Passer `end_depth` en paramètre de ligne de commande
    exported_file = execute_combined_process(end_depth)
    print(f"Process completed, exported file: {exported_file}")
