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
#string

import sys, argparse, random, csv, os, json, string
import numpy as np
from Bio import Phylo
from matplotlib.pyplot import show, savefig
from writeSLiM import writeSLiM
from writeSLiMHPC import writeSLiMHPC

class SLiMTree:

    #Main script to run other commands
    def __init__(self):
        self.read_user_input()
        self.find_fitness_profile()

        clade_data = self.read_clade_data()
        if (self.data_file != None):
            self.data_file.close()

        clade_dict_list = self.read_input_tree(clade_data)

        self.write_slim_code(clade_dict_list)



    #Read input parameters from the user
    def read_user_input(self):

        #Set up starting parameters dictionary
        self.starting_parameters = {}

        #Parse for arguments given by the user
        parser = argparse.ArgumentParser(description='A program to make slim simulations from newick phylogeny files')
        parser.add_argument('-i','--input_tree', nargs = 1, required = True, type = str,
                help = 'tree file in newick format specifying the clades to simulate')
        parser.add_argument('-d','--tree_data_file', nargs = 1, type=argparse.FileType('r'),
                help = 'file specifying population size, mutation rate, etc. for each node, see documentation')
        parser.add_argument('-T', '--tool', type = str, required = True,
                help = 'name of tool you would like to use. Options include SLiM-Tree, SLiM-Tree-HPC. Default = SLiM-Tree')
        parser.add_argument('-p', '--partition', type = str, help = 'partition to run SLiM-Tree HPC on')
        parser.add_argument('-t', '--time', type = str,
                help = 'maximum time to run each simulation for - suggested time is the maximum time available for a partition')


        #Default parameters are somewhat arbitrary and should be changed for each sim
        parser.add_argument('-n','--population_size', help = 'starting population size for the simulation, default = 100', type = int, default = 100)
        parser.add_argument('-v','--mutation_rate', help = 'starting mutation rate for the simulation, default = 2.5e-6', type=float, default = 2.5e-6)
        parser.add_argument('-g','--genome_length', help = 'length of the genome - amino acids, default = 500', type=int, default = 500)
        parser.add_argument('-r','--recombination_rate', help = 'recombination rate, default = 2.5e-8', type=float, default = 2.5e-8)
        parser.add_argument('-b','--burn_in_multiplier', help = 'value to multiply population size by for burn in, default = 10', type=float, default = 10)
        parser.add_argument('-k','--sample_size', help = 'size of sample obtained from each population at output. Input all for whole sample,  default = 10', type=str, default = "10")

        parser.add_argument('-c','--count_subs', type = self.str2bool, default = True, const=True, nargs='?',
                help = 'boolean specifying whether to count substitutions, turning off will speed up sims. default = True')
        parser.add_argument('-o','--output_gens', type = self.str2bool, default = True, const=True, nargs='?',
                help = 'boolean specifying whether to output the generation after every 100th generation. default = True')
        parser.add_argument('-B','--backup', type = self.str2bool, default = True, const=True, nargs='?',
                help = 'boolean specifying whether to backup simulations, turning off will save space. default = True')

        parser.add_argument('-w', '--wright_fisher_model', type = self.str2bool, default=True, const=True, nargs='?',
                help = 'boolean specifying whether this is a wright-fisher model or non-wright-fisher model. default = True')

        parser.add_argument('-G', '--gene_count', type = int, default = 1, help = "Number of genes in the model. Default = 1.")
        parser.add_argument('-C', '--coding_ratio', type = float, default = 1.0, help = "Ratio of the genome which are coding regions as a ratio coding/noncoding. Default = 1.0")

        #Get arguments from user
        arguments = parser.parse_args()

        #Set up tree
        self.input_file = arguments.input_tree[0]

        #Get simulation type and ensure that required arguments are given for simulation type
        self.simulationType = arguments.tool.translate(str.maketrans('', '', string.punctuation)).lower()
        if (self.simulationType == "slimtreehpc" and (arguments.partition == None or arguments.time == None)):
            print("When using SLiM-Tree-HPC, partition and time data must be provided. Closing program.")
            sys.exit(0)

        #Check to make sure gene count and coding ratio are valid
        if (arguments.gene_count < 1 or arguments.gene_count > arguments.genome_length):
            print("Number of genes must be greater than 0 and less than the length of the genome. Closing program.")
            sys.exit(0);

        if (arguments.coding_ratio <= 0 or arguments.coding_ratio > 1.0):
            print("Coding ratio must be greater than 0 and less than or equal to 1. Please re-enter as a ratio (0, 1]. Closing program.")
            sys.exit(0);

        

        #Set up the starting parameters
        self.starting_parameters["mutation_rate"] = arguments.mutation_rate
        self.starting_parameters["population_size"] = arguments.population_size
        self.starting_parameters["genome_length"] = int(arguments.genome_length)
        self.starting_parameters["recombination_rate"] = arguments.recombination_rate
        self.starting_parameters["burn_in"] = arguments.burn_in_multiplier * arguments.population_size
        self.starting_parameters["sample_size"] = arguments.sample_size

        self.starting_parameters["wf_model"] = arguments.wright_fisher_model

        self.starting_parameters["gene_count"] = arguments.gene_count
        self.starting_parameters["coding_ratio"] = arguments.coding_ratio

        self.starting_parameters["partition"] = arguments.partition
        self.starting_parameters["time"] = arguments.time

        self.starting_parameters["count_subs"] = arguments.count_subs
        self.starting_parameters["output_gens"] = arguments.output_gens
        self.starting_parameters["backup"] = arguments.backup

        #Set up the filenames for file io
        input_file_start = os.getcwd() + '/' + self.input_file.split('.')[0]
        self.starting_parameters["tree_filename"] = input_file_start + "_phylogeny.png"
        self.starting_parameters["fasta_filename"] = input_file_start

        if(arguments.tree_data_file != None):
            self.data_file = arguments.tree_data_file[0]
        else:
            self.data_file = None


        #Set up the output of scripts to a single folder
        split_starting_file = input_file_start.split('/')
        output_files_directory = "/".join(split_starting_file[0:(len(split_starting_file)-1)]) + "/slimScripts"
        backup_files_directory = "/".join(split_starting_file[0:(len(split_starting_file)-1)]) + "/backupFiles"


        try:
            os.mkdir(output_files_directory)
            os.mkdir(backup_files_directory)
        except OSError:
            print ("The directory %s already exits, program files will be overwritten" % output_files_directory)

        self.starting_parameters["output_file"] = output_files_directory + "/" + split_starting_file[-1]



        #Save starting parameters and value of theta to a file for later reference
        theta = 4*arguments.mutation_rate*arguments.population_size

        parameter_file = open(input_file_start + "_parameters.txt", 'w')
        parameter_file.write("Simulation parameters\n\n")

        for key, value in self.starting_parameters.items():
            #Don't need to record these filenames as they are not yet complete
            if(key in ['fasta_filename', 'tree_filename', 'output_file']):
                continue

            parameter_file.write('%s:%s\n' % (key, value))

        parameter_file.write("theta: " + str(theta))
        parameter_file.close()


    #Command to take input from user and convert to bool
    #From: https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
    def str2bool(self, v):
        if isinstance(v, bool):
            return v
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')


    #Read fitness profile and stationary distribution data from psi_c50 file, make fitness profiles
    def find_fitness_profile(self):
        #Open stationary and fitness effects csv and make into list
        with open(os.path.dirname(os.path.realpath(__file__)) + '/fitnessDataFiles/table_fitness_profiles.csv', newline='') as fitness_file:
            reader = csv.reader(fitness_file)
            fitness_dist = list(reader)[1:]
            fitness_length = len(fitness_dist)
            #print("Fitness profile length: " + str(fitness_length))
            fitness_profile_nums = random.choices(range(fitness_length),k=self.starting_parameters["genome_length"]-2 )
            #print("Fitness profile nums: " + str(fitness_profile_nums))


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
        self.starting_parameters["fitness_profile_nums"] = fitness_profile_nums
        self.starting_parameters["stationary_distributions"] = stationary_distributions
        self.starting_parameters["fitness_profiles"] = fitness_profiles
        self.starting_parameters["amino_acids"] = amino_acids




    #Read individaul data from the file - to add more parameters modify data_translation_dict
    def read_clade_data(self):

        data_translation_dict = {
            'v': 'mutation_rate',
            'n': 'population_size',
            'r': 'recombination_rate',
            'k': 'sample_size',
            'p': 'partition',
            't': 'time',
            'c': 'count_subs',
            'o': 'output_gens',
            'B': 'backup',
            'mutation_rate': 'mutation_rate',
            'population_size': 'population_size',
            'recombination_rate': 'recombination_rate',
            'sample_size': 'sample_size',
            'partition': 'partition',
            'time': 'time',
            'count_subs': 'count_subs',
            'output_gens': 'output_gens',
            'backup': 'backup'

        }

        if(self.data_file == None):
            return (None)
        else:
            data = {}
            line = self.data_file.readline()

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

                line = self.data_file.readline()
            return(data)



    #Read the phylogenetic tree data given by the user
    def read_input_tree(self, clade_data):

        self.pop_num = 0

        phylogeny = Phylo.read(self.input_file,"newick")

        #Figure out how long the simulation is going to run for
        max_depth = int(max(phylogeny.depths().values())) + self.starting_parameters["burn_in"] + 1
        self.starting_parameters["num_generations"] = max_depth

        #Set up starting parameters and make list of dictionaries of variables
        starting_parameter_dict = {
            "pop_name": None,
            "child_clades" : None,
            "mutation_rate" : self.starting_parameters["mutation_rate"],
            "population_size" : self.starting_parameters["population_size"],
            "recombination_rate" : self.starting_parameters["recombination_rate"],
            "sample_size": self.starting_parameters["sample_size"],
            "partition": self.starting_parameters["partition"],
            "time" : self.starting_parameters["time"],
            "count_subs" : self.starting_parameters["count_subs"],
            "output_gens" : self.starting_parameters["output_gens"],
            "backup" : self.starting_parameters["backup"]
        }

        try:
            clade_dict_list = self.recurse_through_clades(phylogeny.get_nonterminals()[0],
                                                 starting_parameter_dict, clade_data, phylogeny)
        except IndexError:
            print ("Please make sure your input tree is in Newick format. Program closing")
            sys.exit(0)

        #Sort clade dict list by the distance from the start of the simulation so that it works properly
        #in SLiM

        clade_dict_list = sorted(clade_dict_list, key=lambda k: k["dist_from_start"])


        return (clade_dict_list)




    #Recurses through each of the clades in the phylogeny and makes a list of dictionaries of their parameters
    #Recursion is depth first - list is made so that it can be traversed bredth first
    def recurse_through_clades(self, current_clade, parent_clade_dict, clade_data, phylogeny):
        clade_dict_list = self.get_clade_data(current_clade,parent_clade_dict,clade_data, phylogeny)
        clade_dict = clade_dict_list[0]

        #Make list of clades from data
        if (len(clade_dict["child_clades"])==0):
            return (clade_dict_list)
        else:
            #Recurse through all child clades
            list_of_child_clades = []
            for child_clade in clade_dict["child_clades"]:
                child_clade_dict = self.recurse_through_clades(child_clade, clade_dict,clade_data, phylogeny)
                list_of_child_clades = list_of_child_clades + child_clade_dict

            return clade_dict_list + list_of_child_clades




    #Set up data such as clade name, mutation rate, population size, etc. for a clade
    def get_clade_data (self, clade, parent_clade_dict, clade_data, phylogeny):

        #Set up the default parameters based on the parent
        mut_rate = parent_clade_dict["mutation_rate"]
        pop_size = parent_clade_dict["population_size"]
        rec_rate = parent_clade_dict["recombination_rate"]
        samp_size = parent_clade_dict["sample_size"]
        part = parent_clade_dict["partition"]
        time = parent_clade_dict["time"]
        subs = parent_clade_dict["count_subs"]
        gens = parent_clade_dict["output_gens"]
        backup = parent_clade_dict["backup"]

        #Change parameters if specified by user
        if(clade_data != None):
            clade_name = clade.name
            if(clade_name in clade_data.keys()):
                current_clade_data = clade_data[clade_name]

                if('mutation_rate' in current_clade_data.keys()):
                    mut_rate = float(current_clade_data['mutation_rate'])
                if('population_size' in current_clade_data.keys()):
                    pop_size = int(current_clade_data['population_size'])
                if('recombination_rate' in current_clade_data.keys()):
                    rec_rate = float(current_clade_data['recombination_rate'])
                if('sample_size' in current_clade_data.keys()):
                    samp_size = int(current_clade_data['sample_size'])
                if('partition' in current_clade_data.keys()):
                    part = current_clade_data['partition']
                if('time' in current_clade_data.keys()):
                    time = current_clade_data['time']
                if('count_subs' in current_clade_data.keys()):
                    subs = self.str2bool(current_clade_data['count_subs'])
                if('output_gens' in current_clade_data.keys()):
                    gens = self.str2bool(current_clade_data['output_gens'])
                if('backup' in current_clade_data.keys()):
                    backup = self.str2bool(current_clade_data['backup'])


        #Figure out what population name is for self and assign clade name appropriately
        self.pop_num += 1
        pop_name = "p" + str(self.pop_num)

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
            "end_dist": self.starting_parameters["burn_in"]  + phylogeny.distance(clade),
            "terminal_clade" : clade.clades == [],
            "last_child_clade" : last_child_clade,
            "sample_size": samp_size,
            "partition": part,
            "time" : time,
            "count_subs" : subs,
            "output_gens" : gens,
            "backup" : backup
        }

        return [clade_dict]


    #Using the dictionary of the data from the clades, write slim code to represent the phylogeny
    def write_slim_code (self, clade_dict_list):

        #Open SLiM writer based on tool type and write the initialize statement
        print(self.simulationType)

        if(self.simulationType == "slimtree"):
            SLiM_Writer = writeSLiM(self.starting_parameters)
        elif(self.simulationType == "slimtreehpc"):
            SLiM_Writer = writeSLiMHPC(self.starting_parameters)
        else:
            print ("Invalid tool type. Please specify a tool as SLiM-Tree or SLiM-Tree-HPC. Program closing")
            sys.exit(0)

        #Write a script for each clade which will be run sequentially
        if (self.starting_parameters["wf_model"]): #If this is a Wright-Fisher model, use a different write_subpop function
            for clade in clade_dict_list:
                SLiM_Writer.write_subpop(clade)
        else:
            for clade in clade_dict_list:
                SLiM_Writer.write_subpop_nonwf(clade);

        #Start the SLiM code to run
        if(self.simulationType == "slimtree"):
            SLiM_Writer.close_file()
            os.system("slim \"" + self.starting_parameters["output_file"] + "_p1.slim\"")
        elif(self.simulationType == "slimtreehpc"):
            os.system("sbatch \"" + self.starting_parameters["output_file"] + "_p1.sh\"")



if __name__=='__main__':
    SLiMTree()
