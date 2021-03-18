#Program to take in a tree in Newick file and output SLiM code to run a simulation from that tree

#Required packages:
#sys
#argparse
#BioPython
#matplotlib
#random
#csv
#numpy
#os
#json

import sys, argparse, random, csv, os, json
import numpy as np
from Bio import Phylo
from matplotlib.pyplot import show, savefig
from writeSLiM import writeSLiM


#Read input parameters from the user
def read_user_input():
    global input_file
    global data_file

    global starting_parameters


    #Set up starting parameters dictionary 
    starting_parameters = {}

    
    parser = argparse.ArgumentParser(description='A program to make slim simulations from newick phylogeny files')
    parser.add_argument('-i','--input_tree', nargs = 1, required = True, type = str,
            help = 'tree file in newick format specifying the clades to simulate')
    parser.add_argument('-d','--tree_data_file', nargs = 1, type=argparse.FileType('r'),
            help = 'file specifying population size, mutation rate, etc. for each node, see documentation')
    
    #Default parameters give a default theta of 0.01
    parser.add_argument('-n','--population_size', help = 'starting population size for the simulation, default = 500', type = int, default = 500)
    parser.add_argument('-v','--mutation_rate', help = 'starting mutation rate for the simulation, default = 0.000001', type=float, default = 0.000001)
    parser.add_argument('-g','--genome_length', help = 'length of the genome - amino acids, default = 3.33e5', type=int, default = 3.33e5)
    parser.add_argument('-r','--recombination_rate', help = 'recombination rate, default = 1e-8', type=float, default = 1e-8)
    parser.add_argument('-b','--burn_in_multiplier', help = 'value to multiply popsize by for burn in, default = 10', type=float, default = 10)
    parser.add_argument('-k','--sample_size', help = 'size of sample obtained from each population at output. Input all for whole sample,  default = 10', type=str, default = "10")

    #Set up important variables
    arguments = parser.parse_args()

    input_file = arguments.input_tree[0]

    #Set up the starting parameters
    starting_parameters["mutation_rate"] = arguments.mutation_rate
    starting_parameters["population_size"] = arguments.population_size
    starting_parameters["genome_length"] = int(arguments.genome_length)
    starting_parameters["recombination_rate"] = arguments.recombination_rate
    starting_parameters["burn_in"] = arguments.burn_in_multiplier * arguments.population_size
    starting_parameters["sample_size"] = arguments.sample_size

    #Set up the filenames for file io
    input_file_start = os.getcwd() + '/' + input_file.split('.')[0] 
    starting_parameters["tree_filename"] = input_file_start + "_phylogeny.png"
    starting_parameters["fasta_filename"] = input_file_start

    if(arguments.tree_data_file != None):
        data_file = arguments.tree_data_file[0]
    else:
        data_file = None  


    #Set up the output of scripts to a single folder
    split_starting_file = input_file_start.split('/')
    output_files_directory = "/".join(split_starting_file[0:(len(split_starting_file)-1)]) + "/slimScripts/"
    backup_files_directory = "/".join(split_starting_file[0:(len(split_starting_file)-1)]) + "/backupFiles/"
     

    try:
        os.mkdir(output_files_directory)
        os.mkdir(backup_files_directory)
    except OSError:
        print ("The directory %s already exits, program files will be overwritten" % output_files_directory)

    starting_parameters["output_file"] = output_files_directory + "/" + split_starting_file[-1]


 
    #Save starting parameters and value of theta to a file for later reference
    theta = 4*arguments.mutation_rate*arguments.population_size

    parameter_file = open(input_file_start + "_parameters.txt", 'w')
    parameter_file.write("Simulation parameters\n\n")
    
    for key, value in starting_parameters.items():
        #Don't need to record these filenames as they are not yet complete
        if(key in ['fasta_filename', 'tree_filename', 'output_file']):
            continue
        
        parameter_file.write('%s:%s\n' % (key, value))
        
    parameter_file.write("theta: " + str(theta))
    parameter_file.close()




#Read fitness profile and stationary distribution data from psi_c50 file, make fitness profiles
def find_fitness_profile():
    global starting_parameters
    
    fitness_profile_nums = random.choices(range(50),k=starting_parameters["genome_length"]-2 )
    #print("Fitness profile nums:" + str(fitness_profile_nums))

    #Open stationary and fitness effects csv and make into list
    with open(os.path.dirname(os.path.realpath(__file__)) + '/fitnessDataFiles/table_fitness_profiles.csv', newline='') as fitness_file:
        reader = csv.reader(fitness_file)
        fitness_dist = list(reader)[1:]
    

    with open(os.path.dirname(os.path.realpath(__file__)) + '/fitnessDataFiles/table_stationary_distributions.csv', newline='') as stationary_file:
        reader = csv.reader(stationary_file)
        stationary_dist = list(reader)[1:]


    #Find stationary distribution of each of the 53 fitness profiles
    stationary_distributions = np.array(stationary_dist)
    stationary_distributions = stationary_distributions.astype(float)


    #Find fitness effects for each amino acid
    fitness_dist = np.array(fitness_dist).astype(float)
    fitness_dist = fitness_dist.tolist()
    amino_acids = (["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N"] +
                           ["P", "Q", "R", "S", "T", "V", "W", "Y", "X"])
    fitness_profiles = {}
    
    for amino_acid_num in range(len(amino_acids)):
            aa = amino_acids[amino_acid_num]
            fitness_profiles[aa] = fitness_dist[amino_acid_num]


    #Set up distributions in starting parameters
    starting_parameters["fitness_profile_nums"] = fitness_profile_nums
    starting_parameters["stationary_distributions"] = stationary_distributions
    starting_parameters["fitness_profiles"] = fitness_profiles
    starting_parameters["amino_acids"] = amino_acids




#Read individaul data from the file - to add more parameters modify data_translation_dict
def read_clade_data():
    global data_file

    data_translation_dict = {
        'v': 'mutation_rate',
        'n': 'population_size',
        'r': 'recombination_rate',
        'mutation_rate': 'mutation_rate',
        'population_size': 'population_size',
        'recombination_rate': 'recombination_rate'
        
    }

    if(data_file == None):
        return (None)
    else:
        data = {}
        line = data_file.readline()

        while (line != ''):
            line = line.split('\n')[0]
            if (line == ''):
                pass
            elif(line[0] == '@'):
                data_for_node = {}
                data[line[1:]] = data_for_node
            elif(line[0] == '-'):
                data_label = line[1]
                data_value = line.split(' ')[1]
                
                data_for_node[data_translation_dict[data_label]] = data_value
            
            line = data_file.readline()
        return(data)



#Read the phylogenetic tree data given by the user
def read_input_tree(clade_data):
    global input_file
    global starting_parameters
    global pop_num

    pop_num = 0
    
    phylogeny = Phylo.read(input_file,"newick")

    #Figure out how long the simulation is going to run for
    max_depth = int(max(phylogeny.depths().values())) + starting_parameters["burn_in"] + 1
    starting_parameters["num_generations"] = max_depth

    #Set up starting parameters and make list of dictionaries of variables
    starting_parameter_dict = {
        "pop_name": None,
        "child_clades" : None,
        "mutation_rate" : starting_parameters["mutation_rate"],
        "population_size" : starting_parameters["population_size"],
        "recombination_rate" : starting_parameters["recombination_rate"]
    }

    try:
        clade_dict_list = recurse_through_clades(phylogeny.get_nonterminals()[0],
                                             starting_parameter_dict, clade_data, phylogeny)
    except IndexError:
        print ("Please make sure your input tree is in Newick format. Program closing")
        sys.exit(0)

    #Sort clade dict list by the distance from the start of the simulation so that it works properly
    #in SLiM

    clade_dict_list = sorted(clade_dict_list, key=lambda k: k["dist_from_start"]) 
    
##    print(clade_dict_list)

    
    return (clade_dict_list)




#Recurses through each of the clades in the phylogeny and makes a list of dictionaries of their parameters
#Recursion is depth first - list is made so that it can be traversed bredth first
def recurse_through_clades(current_clade, parent_clade_dict, clade_data, phylogeny):
    clade_dict_list = get_clade_data(current_clade,parent_clade_dict,clade_data, phylogeny)
    clade_dict = clade_dict_list[0]
    
    #Make list of clades from data
    if (len(clade_dict["child_clades"])==0):
        return (clade_dict_list)
    else:
        #Recurse through all child clades
        list_of_child_clades = []
        for child_clade in clade_dict["child_clades"]:
            child_clade_dict = recurse_through_clades(child_clade, clade_dict,clade_data, phylogeny)
            list_of_child_clades = list_of_child_clades + child_clade_dict
        
        return clade_dict_list + list_of_child_clades       




#Set up data such as clade name, mutation rate, population size, etc. for a clade
def get_clade_data (clade, parent_clade_dict, clade_data, phylogeny):
    global pop_num
    global starting_parameters

    #Set up the default mutation rate, population size and recombination rate based on the parent
    mut_rate = parent_clade_dict["mutation_rate"]
    pop_size = parent_clade_dict["population_size"]
    rec_rate = parent_clade_dict["recombination_rate"]

    #Change mutation rate, population size and recombination rate if specified by user
    if(clade_data != None):
        clade_name = clade.name
        if(clade_name in clade_data.keys()):
            current_clade_data = clade_data[clade_name]

            if('mutation_rate' in current_clade_data.keys()):
                mut_rate = float(current_clade_data['mutation_rate'])
            if('population_size' in current_clade_data.keys()):
                pop_size = int(current_clade_data['population_size'])
            if('recombination_rate' in current_clade_data.keys()):
                rec_rate = int(current_clade_data['recombination_rate'])

    #Figure out what population name is for self and assign clade name appropriately
    pop_num += 1
    pop_name = "p" + str(pop_num)

    if(clade.name != None):
        clade.name = pop_name + ": "+ clade.name
    else:
        clade.name = pop_name

    #Figure out when the population needs to be formed
    dist_from_parent = clade.branch_length
    if(dist_from_parent == None):
        dist_from_start = 0
    else:
        dist_from_start = parent_clade_dict["end_dist"]


    #Determine whether population belongs to the last child clade - allows for removal of extraneous data
    parents_children = parent_clade_dict["child_clades"]

    if(parents_children == None):
        last_child_clade = False
    else:
        last_child_clade = clade == parents_children[-1]
    
    

    #Set up the dictionary of values for the clade
    clade_dict = {
        "pop_name": pop_name,
        "parent_pop_name" : parent_clade_dict["pop_name"],
        "child_clades" : clade.clades, 
        "mutation_rate" : mut_rate,
        "population_size" : pop_size,
        "recombination_rate": rec_rate, 
        "dist_from_start" : dist_from_start,
        "end_dist": starting_parameters["burn_in"]  + phylogeny.distance(clade),
        "terminal_clade" : clade.clades == [],
        "last_child_clade" : last_child_clade
    }
    
    return [clade_dict]


#Using the dictionary of the data from the clades, write slim code to represent the phylogeny
def write_slim_code (clade_dict_list):
    global starting_parameters

    #Open SLiM writer and write the initialize statement
    SLiM_Writer = writeSLiM(starting_parameters)

    #Write a script for each clade which will be run sequentially
    for clade in clade_dict_list:
        SLiM_Writer.write_subpop(clade)

    #Start the SLiM code to run
    os.system("sbatch \"" + starting_parameters["output_file"] + "_p1.sh\"")





#Main script to run other commands
def main():
    read_user_input()
    find_fitness_profile()

    global data_file
    clade_data = read_clade_data()
    if (data_file != None):
        data_file.close()
    
    clade_dict_list = read_input_tree(clade_data)

    write_slim_code(clade_dict_list)
    
    




if __name__=='__main__':
    main()
