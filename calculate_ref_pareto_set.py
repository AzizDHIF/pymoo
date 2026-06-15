import numpy as np

from pymoo.config import Config
Config.warnings['not_compiled'] = False

from pymoo.core.problem import ElementwiseProblem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize

from pymoo.operators.sampling.rnd import BinaryRandomSampling
from pymoo.operators.crossover.pntx import TwoPointCrossover
from pymoo.operators.mutation.bitflip import BitflipMutation

import sys
import subprocess
import os 
WORK_DIR = os.path.dirname(os.path.abspath(__file__))

def make_negative_pareto(filename):
     with open(filename, "r") as  f:
          lines = f.readlines()
    
     for i  in range(len(lines)):
        line=lines[i]
        line=line.strip().split()
        
        line[0] = str(-float(line[0]))
        line[1] = str(-float(line[1]))
        line[2] = str(-float(line[2]))

        lines[i]=" ".join(line) + "\n"
      
     with open(filename, "w") as r:
          r.writelines(lines)
          
          
def read_instance(filename):

    with open(filename, "r") as f:
        lines = [l.strip() for l in f if l.strip()]

    n_obj, n_items = map(int, lines[0].split())

    capacities = np.zeros(n_obj)
    weights = np.zeros((n_obj, n_items))
    profits = np.zeros((n_obj, n_items))

    idx = 1

    for obj in range(n_obj):

        capacities[obj] = float(lines[idx])
        idx += 1

        for item in range(n_items):

            idx += 1  # item_x

            weights[obj, item] = float(lines[idx])
            idx += 1

            profits[obj, item] = float(lines[idx])
            idx += 1

    return n_items, n_obj, capacities, weights, profits

def run_aco(executable, args=[]):
    run_cmd = [f"./{executable}.exe"] + args
    print(f"Exécution : {' '.join(run_cmd)}")
    
    result = subprocess.run(run_cmd, cwd=WORK_DIR)  
    
    if result.returncode != 0:
        print(f"Erreur à l'exécution (code {result.returncode})")
        sys.exit(1)

def extraire_pareto_sets(fichier_entree, dossier_sortie):
    import os

    os.makedirs(dossier_sortie, exist_ok=True)

    fichier_sortie = None
    compteur = 0

    with open(fichier_entree, 'r') as f:
        for ligne in f:
            ligne = ligne.strip()
            if not ligne:
                continue
            if ligne.startswith("mcycle"):
                if fichier_sortie:
                    fichier_sortie.close()
                compteur += 1
                chemin = os.path.join(dossier_sortie, f"pareto_front_aco_seed_{compteur}.txt")
                fichier_sortie = open(chemin, 'w')
                continue

            # Écrire uniquement les lignes de points (3 flottants)
            parties = ligne.split()
            if len(parties) == 3:
                try:
                    parties_negatives=[-1*float(x) for x in parties]
                    ligne_negative = "\t".join(str(v) for v in parties_negatives)

                    fichier_sortie.write(ligne_negative + "\n")
                except ValueError:
                    pass

    if fichier_sortie:
        fichier_sortie.close()

    print(f"{compteur} fichiers créés dans '{dossier_sortie}'")
class MOMDKP(ElementwiseProblem):

    def __init__(self, filename):

        self.n_items, self.n_obj, self.capacities, \
        self.weights, self.profits = read_instance(filename)

        super().__init__(
            n_var=self.n_items,
            n_obj=self.n_obj,
            n_ieq_constr=self.n_obj,
            xl=0,
            xu=1
        )

    def _evaluate(self, x, out, *args, **kwargs):

        F = []

        for k in range(self.n_obj):
            F.append(-np.sum(self.profits[k] * x))

        G = []

        for k in range(self.n_obj):
            G.append(
                np.sum(self.weights[k] * x)
                - self.capacities[k]
            )

        out["F"] = np.array(F)
        out["G"] = np.array(G)



from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.util.ref_dirs import get_reference_directions




from pymoo.algorithms.moo.moead import MOEAD
from pymoo.util.ref_dirs import get_reference_directions

ref_dirs = get_reference_directions(
    "das-dennis",
    3,
    n_partitions=12
)




algorithms = {
    "NSGA2": lambda: NSGA2(
        pop_size=100,
        sampling=BinaryRandomSampling(),
        crossover=TwoPointCrossover(),
        mutation=BitflipMutation()
    ),

    "NSGA3": lambda: NSGA3(
        pop_size=100,
        ref_dirs=ref_dirs,
        sampling=BinaryRandomSampling(),
        crossover=TwoPointCrossover(),
        mutation=BitflipMutation()
    ),



}
def append_pareto_to_res(pareto_file, res_file):
    with open(pareto_file, 'r') as f:
        lines = f.readlines()
    with open(res_file, 'a') as f:
        f.writelines(lines)

if __name__ == "__main__":
    mood=sys.argv[1]
    if mood=="train":
        datasets=["pymoo\\dataset\\mood_train_dataset\\dataset_0_instance_100_items_3_objectifs.txt","pymoo\\dataset\\mood_train_dataset\\dataset_1_instance_100_items_3_objectifs.txt","pymoo\\dataset\\mood_train_dataset\\dataset_2_instance_100_items_3_objectifs.txt","pymoo\\dataset\\mood_train_dataset\\dataset_3_instance_100_items_3_objectifs.txt","pymoo\\dataset\\mood_train_dataset\\dataset_4_instance_100_items_3_objectifs.txt"]

        print("Exécution des algorithmes NSGA2 et NSGA3 et WeightACO sur les instances du dataset...")
        for algo_name, algo_factory in algorithms.items():
            
            for i in range(5):
                    for seed in [1, 2, 3, 4, 5,6,7,8,9,10]:
                        algo = algo_factory()

                        print(f"\n=== {algo.__class__.__name__} | instance {i} ==")

                        problem = MOMDKP(f"pymoo\\dataset\\mood_train_dataset\\dataset_{i}_instance_100_items_3_objectifs.txt")


                        res = minimize(
                            problem,
                            algo,
                            ('n_gen', 300),
                            seed=seed,
                            verbose=True
                        )

                        np.savetxt(rf"pymoo\\dataset\\mood_train_dataset\\ref_dataset_{i}\\pareto_front_{algo.__class__.__name__.lower()}_seed_{seed}.txt", res.F, fmt="%.6f")


        print("Lancement pour l'ACO ...")

        for dataset in datasets:
                    run_aco("WeightACO_100items",args=[dataset])


        #extraire les sets de pareto à partir des résultats de l'ACO

        aco_results=[("results_train_dataset_0.txt", "pymoo\\dataset\\mood_train_dataset\\ref_dataset_0"),
                        ("results_train_dataset_1.txt", "pymoo\\dataset\\mood_train_dataset\\ref_dataset_1"),
                        ("results_train_dataset_2.txt", "pymoo\\dataset\\mood_train_dataset\\ref_dataset_2"),
                        ("results_train_dataset_3.txt", "pymoo\\dataset\\mood_train_dataset\\ref_dataset_3"),
                        ("results_train_dataset_4.txt", "pymoo\\dataset\\mood_train_dataset\\ref_dataset_4")]
            
        for (algo_result_file, pareto_output_dir) in aco_results:
            extraire_pareto_sets(algo_result_file, pareto_output_dir)
        

        #rendre les front de pareto produits par ACO négative

        for seed in [1, 2, 3, 4, 5,6,7,8,9,10]:
           for i in range(5):
                                
             pareto_file = f"pymoo\\dataset\\mood_train_dataset\\ref_dataset_{i}\\pareto_front_aco_seed_{seed}.txt"

             make_negative_pareto(pareto_file)


        print("\nTous les algorithmes ont été exécutés avec succès. concaténation des front de pareto...")


        for i in [0, 1, 2, 3, 4]:
                
                    for algo in ["nsga2", "nsga3","aco"]:
                        for seed in [1, 2, 3, 4, 5,6,7,8,9,10]:
                                
                                pareto_file = f"pymoo\\dataset\\mood_train_dataset\\ref_dataset_{i}\\pareto_front_{algo}_seed_{seed}.txt"
                                res_file = f"pymoo\\dataset\\mood_train_dataset\\ref_dataset_{i}\\final_pareto_file.txt"
                                append_pareto_to_res(pareto_file, res_file)

        



        print("\nTous les front de pareto ont été concaténés. Exécution de l'outil nondominated.exe pour filtrer les solutions non dominées...")


        for i in [0, 1, 2, 3, 4]:
            
                res_file = f"pymoo\\dataset\\mood_train_dataset\\ref_dataset_{i}\\final_pareto_file.txt"
                run_cmd = [
                "./nondominated.exe",
                "--union",
                "--filter",
                res_file]
                print(f"Exécution : {' '.join(run_cmd)}")
                
                result = subprocess.run(run_cmd, cwd=WORK_DIR)  
                
                if result.returncode != 0:
                    print(f"Erreur à l'exécution (code {result.returncode})")
                    sys.exit(1)

        print("\nToutes les solutions non dominées ont été filtrées avec succès. Vérification  de la non-dominance des solutions filtrées...")


        for i in [0, 1, 2, 3, 4]:
                res_file = f"pymoo\\dataset\\mood_train_dataset\\ref_dataset_{i}\\final_pareto_file.txt_dat"
                run_cmd = [
                "./nondominated.exe",
                "--verbose",
                res_file]
                print(f"Exécution : {' '.join(run_cmd)}")
                
                result = subprocess.run(run_cmd, cwd=WORK_DIR)  
                
                if result.returncode != 0:
                    print(f"Erreur à l'exécution (code {result.returncode})")
                    sys.exit(1)
    
    else:
        from itertools import product

        datasets=[f"pymoo\\dataset\\mood_val_dataset\\dataset_{i}_instance_{nb_items}_items_3_objectifs.txt" for i,nb_items in product([0,1,2,3,4],[100,300,500])]
        
        print("Exécution des algorithmes NSGA2 et NSGA3 et WeightACO sur les instances du dataset...")
        for algo_name, algo_factory in algorithms.items():
            
            for i in range(5):
                    for seed in [1, 2, 3, 4, 5,6,7,8,9,10]:
                        algo = algo_factory()

                        print(f"\n=== {algo.__class__.__name__} | instance {i} ==")

                        problem = MOMDKP(f"pymoo\\dataset\\mood_train_dataset\\dataset_{i}_instance_100_items_3_objectifs.txt")


                        res = minimize(
                            problem,
                            algo,
                            ('n_gen', 300),
                            seed=seed,
                            verbose=True
                        )

                        np.savetxt(rf"pymoo\\dataset\\mood_train_dataset\\ref_dataset_{i}\\pareto_front_{algo.__class__.__name__.lower()}_seed_{seed}.txt", res.F, fmt="%.6f")


        print("Lancement pour l'ACO ...")

        for dataset in datasets:
                    run_aco("WeightACO_100items",args=[dataset])


        #extraire les sets de pareto à partir des résultats de l'ACO

        aco_results=[("results_train_dataset_0.txt", "pymoo\\dataset\\mood_train_dataset\\ref_dataset_0"),
                        ("results_train_dataset_1.txt", "pymoo\\dataset\\mood_train_dataset\\ref_dataset_1"),
                        ("results_train_dataset_2.txt", "pymoo\\dataset\\mood_train_dataset\\ref_dataset_2"),
                        ("results_train_dataset_3.txt", "pymoo\\dataset\\mood_train_dataset\\ref_dataset_3"),
                        ("results_train_dataset_4.txt", "pymoo\\dataset\\mood_train_dataset\\ref_dataset_4")]
            
        for (algo_result_file, pareto_output_dir) in aco_results:
            extraire_pareto_sets(algo_result_file, pareto_output_dir)
        

        #rendre les front de pareto produits par ACO négative

        for seed in [1, 2, 3, 4, 5,6,7,8,9,10]:
           for i in range(5):
                                
             pareto_file = f"pymoo\\dataset\\mood_train_dataset\\ref_dataset_{i}\\pareto_front_aco_seed_{seed}.txt"

             make_negative_pareto(pareto_file)


        print("\nTous les algorithmes ont été exécutés avec succès. concaténation des front de pareto...")


        for i in [0, 1, 2, 3, 4]:
                
                    for algo in ["nsga2", "nsga3","aco"]:
                        for seed in [1, 2, 3, 4, 5,6,7,8,9,10]:
                                
                                pareto_file = f"pymoo\\dataset\\mood_train_dataset\\ref_dataset_{i}\\pareto_front_{algo}_seed_{seed}.txt"
                                res_file = f"pymoo\\dataset\\mood_train_dataset\\ref_dataset_{i}\\final_pareto_file.txt"
                                append_pareto_to_res(pareto_file, res_file)

        



        print("\nTous les front de pareto ont été concaténés. Exécution de l'outil nondominated.exe pour filtrer les solutions non dominées...")


        for i in [0, 1, 2, 3, 4]:
            
                res_file = f"pymoo\\dataset\\mood_train_dataset\\ref_dataset_{i}\\final_pareto_file.txt"
                run_cmd = [
                "./nondominated.exe",
                "--union",
                "--filter",
                res_file]
                print(f"Exécution : {' '.join(run_cmd)}")
                
                result = subprocess.run(run_cmd, cwd=WORK_DIR)  
                
                if result.returncode != 0:
                    print(f"Erreur à l'exécution (code {result.returncode})")
                    sys.exit(1)

        print("\nToutes les solutions non dominées ont été filtrées avec succès. Vérification  de la non-dominance des solutions filtrées...")


        for i in [0, 1, 2, 3, 4]:
                res_file = f"pymoo\\dataset\\mood_train_dataset\\ref_dataset_{i}\\final_pareto_file.txt_dat"
                run_cmd = [
                "./nondominated.exe",
                "--verbose",
                res_file]
                print(f"Exécution : {' '.join(run_cmd)}")
                
                result = subprocess.run(run_cmd, cwd=WORK_DIR)  
                
                if result.returncode != 0:
                    print(f"Erreur à l'exécution (code {result.returncode})")
                    sys.exit(1)"""
         
         

