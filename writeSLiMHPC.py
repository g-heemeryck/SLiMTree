#Program to write SLiMTree code for a high performance computing cluster


import random, csv, os
from writeSLiM import writeSLiM



class writeSLiMHPC(writeSLiM):
    #Write code to add a new population to the simulation by writing a new script for that population
    def write_subpop(self, population_parameters):
        pop_name = population_parameters["pop_name"]


        #Create a new script and batch file for the population
        batch_file = open(self.general_output_filename + "_" + pop_name + ".sh", "w")
        batch_file.write("#!/bin/sh\n\n#SBATCH -J SLiM_Simulation_" + pop_name + "\n#SBATCH -t " + population_parameters["time"] +
                "\n#SBATCH -p "  + population_parameters["partition"] + "\n#SBATCH -o " + pop_name + ".out\n#SBATCH -e " + pop_name +
                ".err\n#SBATCH -n 1" + "\n\nslim " + self.general_output_filename + "_" + pop_name +".slim")
        batch_file.close()

        self.output_file = open(self.general_output_filename + "_" + pop_name + ".slim" , "w")

        #Set up the initialize and fitness functions for the new script
        super().write_initialize(population_parameters)
        super().write_fitness()

        #Write the commands that are run for every simulation and the starting population
        self.write_repeated_commands(population_parameters)
        self.write_start_pop(population_parameters)


        #Finish writing the script
        self.write_end_sim(population_parameters)
        self.output_file.close()




    def write_subpop_nonwf(self, population_parameters):
        pop_name = population_parameters["pop_name"]
    
        #Create a new script and batch file for the population
        batch_file = open(self.general_output_filename + "_" + pop_name + ".sh", "w")
        batch_file.write("#!/bin/sh\n\n#SBATCH -J SLiM_Simulation_" + pop_name + "\n#SBATCH -t " + population_parameters["time"] +
                "\n#SBATCH -p "  + population_parameters["partition"] + "\n#SBATCH -o " + pop_name + ".out\n#SBATCH -e " + pop_name
                +".err\n#SBATCH -n 1" + "\n\nslim " + self.general_output_filename + "_" + pop_name+".slim")
        batch_file.close()

        self.output_file = open(self.general_output_filename + "_" + pop_name + ".slim" , "w")

        #Set up the initialize and fitness functions for the new script
        super().write_initialize(population_parameters)
        super().write_fitness()

        super().write_reproduction()

        #Write the commands that are run for every simulation and the starting population
        self.write_repeated_commands(population_parameters)
        self.write_start_pop(population_parameters)

        self.write_early_function(int(population_parameters["dist_from_start"]) +1, int(population_parameters["end_dist"]), population_parameters)


        #Finish writing the script
        self.write_end_sim(population_parameters)
        self.output_file.close()


    #Write code to set up the starting population for each simulation. If first population, population established, otherwise starting population is loaded
    def write_start_pop(self, population_parameters):

        pop_string = (str(int(population_parameters["dist_from_start"])+1) + " late() {" +
                    "\n\tsetup_fitness();")

        #If first population make the population, otherwise load from the parent
        if(population_parameters["parent_pop_name"] == None):
            pop_string += ("\n\twriteFile(\"" + self.fasta_filename + "_aa.fasta\", \"\", append = F);" +
                                   "\n\twriteFile(\"" + self.fasta_filename + "_nuc.fasta\", \"\", append = F);" +
                                   "\n\tsim.addSubpop(\"p1\", " + str(population_parameters["population_size"]) + ");" )

            #Write code to start a fixed state from the starting nucleotide sequence
            pop_string += "\n\tsim.setValue(\"fixations\", strsplit(sim.chromosome.ancestralNucleotides(),sep = \"\"));"
        else:

            #Set appropriate starting population size
            if (self.model_type == True):
                pop_string += ("\n\tsim.readFromPopulationFile(\"" + population_parameters["parent_pop_name"]  + ".txt\");")
                pop_string += ("\n\tp1.setSubpopulationSize(" + str(population_parameters["population_size"]) + ");")
            else:
                #If a non-WF model, take half of the individuals from the parent population to represent the population split according to tags assigned in previous generation.
                if (population_parameters["last_child_clade"]):
                    #Have population tag 1 have fitness 0.0 so they won't influence next generation
                    pop_string += ("\n\tsim.readFromPopulationFile(\"" + population_parameters["parent_pop_name"]  + "_2.txt\");")
                    pop_string += ("\n\tp2.removeSubpopulation();")
                else:
                    #Have population tag 2 have fitness 0.0 so they won't influence next generation.
                    pop_string += ("\n\tsim.readFromPopulationFile(\"" + population_parameters["parent_pop_name"]  + "_1.txt\");")
                    pop_string += ("\n\tp2.removeSubpopulation();")

            #Load population into the end of the parent population's script to start this script when parent's finishes
            parent_output_file = open(self.general_output_filename + "_" + population_parameters["parent_pop_name"] + ".slim" , "a")
            parent_output_file.write("\n\tsystem(\"sbatch \\\"" + self.general_output_filename + "_" + population_parameters["pop_name"] + ".sh\\\"\");")

            if(population_parameters["last_child_clade"]):
                parent_output_file.write("\n}")
            
            parent_output_file.close()

            #Write code to import in the prevouisly fixed state
            pop_string += ("sim.setValue(\"fixations\", strsplit(readFile(\""+ population_parameters["parent_pop_name"] +
                           "_fixed_mutations.txt\"), sep = \"\"));")

        #At the start of the sim there are no fixations counted
        pop_string += "\n\tsim.setValue(\"fixations_counted\", 0);"
        pop_string += "\n}\n\n\n"

        self.output_file.write(pop_string)


    #Write code for early functions in nonWF models.
    def write_early_function(self, start_dist, end_dist, population_parameters):
            early_event = str(int(population_parameters["dist_from_start"]) + 2) + ":" + str(int(population_parameters["end_dist"]) + 1) + " early(){"
            early_event += "\n\tp1.fitnessScaling = " + str(int(population_parameters["population_size"])) + "/ p1.individualCount;" + "\n}\n\n\n"
            self.output_file.write(early_event)


    #Write code to count substitutions, make a backup and count generations
    def write_repeated_commands(self, population_parameters):
        #Set up variables for repeated commands
        start_dist = int(population_parameters["dist_from_start"])+1
        end_dist = int(population_parameters["end_dist"])
        pop_name =  population_parameters["pop_name"]
        
        repeated_commands_string = str(start_dist) +":" + str(end_dist) + "late () {"

        #Write a command to count the substitutions (identity by state)
        if (population_parameters["count_subs"]):
            repeated_commands_string += ("\n\tif(length(sim.mutations)!= 0){"
                        "\n\t\tancestral_genome = sim.getValue(\"fixations_p1\");" +
                        "\n\t\tcompare_genome = strsplit(p1.genomes[0].nucleotides(), sep = \'\');"+
                        "\n\t\tfixed_nucs = rep(T, length(compare_genome));" +
                        "\n\n\t\tfor (genome in (p1.genomes)){" +
                        "\n\t\t\tsame_nucs = (compare_genome == strsplit(genome.nucleotides(), sep = \'\'));" +
                        "\n\t\t\tfixed_nucs = (fixed_nucs & same_nucs);\n\t\t}" +
                        "\n\n\t\tdifferent_muts = (ancestral_genome != compare_genome);" +
                        "\n\t\tnew_fixations = different_muts & fixed_nucs;" +
                        "\n\t\tsim.setValue(\"fixations_counted_p1\", sim.getValue(\"fixations_counted_p1\") + sum(new_fixations));" +
                        "\n\n\t\tancestral_genome[new_fixations] = compare_genome[new_fixations];" +
                        "\n\t\tsim.setValue(\"fixations_p1\", ancestral_genome);\n\t};")

        #Write a command to output when every 100th generation has passed
        if(population_parameters["output_gens"]):
            repeated_commands_string += "\n\n\tif (sim.generation%100 == 0) {\n\t\tcatn(sim.generation);\n\t};"


        #Write a command to write a backup of all individuals after every 100 generations
        if (population_parameters["backup"]):
             repeated_commands_string += ("\n\n\tif (sim.generation%100 == 0) {" +
                        "\n\t\twriteFile(\"" + os.getcwd()+ "/backupFiles/" + pop_name + ".fasta\"," +
                        "(\">parent_ancestral_to_load\\n\" + sim.chromosome.ancestralNucleotides()));" +
                        "\n\t\tsim.outputFull(\"" + os.getcwd()+ "/backupFiles/" + pop_name + ".txt\");\n\t};")

        repeated_commands_string += "\n}\n\n\n"

        self.output_file.write(repeated_commands_string)


    #Write the closing statements to end the simulation and either allow for starting of subsequent simulations or output data (depending on terminal status)
    def write_end_sim(self, population_parameters):

        end_population_string = str(int(population_parameters["end_dist"]) + 1) + " late() {"

        #If terminal clade output data otherwise create data to be loaded into scripts of the clades children
        if(population_parameters["terminal_clade"]):
            end_population_string += super().write_terminal_output(population_parameters)
        else:
            if (self.model_type == False):
                #Tag each individual with either 1 or 2 to go into different subpopulations. Should be split evenly.
                end_population_string += "\n\tp1.individuals.tag = 0;\n\tsample(p1.individuals, integerDiv(p1.individualCount, 2)).tag = 1;\n\tp1.individuals[p1.individuals.tag == 0].tag = 2;\n\tsim.addSubpop(\"p2\", 0);"
                end_population_string += "\n\tp2.takeMigrants(p1.individuals[p1.individuals.tag == 2]);\n\tsim.outputFull(\""+ population_parameters["pop_name"] +"_1.txt\");\n\tp1.takeMigrants(p2.individuals);\n\tp2.takeMigrants(p1.individuals[p1.individuals.tag == 1]);"
                end_population_string += "\n\tsim.outputFull(\""+ population_parameters["pop_name"] +"_2.txt\");\n\tp1.takeMigrants(p2.individuals);\n\tp2.removeSubpopulation();"
            else:
                end_population_string += "\n\tsim.outputFull(\"" + population_parameters["pop_name"] + ".txt\");"

            end_population_string += ("\n\twriteFile(\"" + population_parameters["pop_name"] +
                                      ".fasta\", (\">parent_ancestral_to_load\\n\" + sim.chromosome.ancestralNucleotides()));")

        #Scripting to end the simulation and write the fixed mutations
        end_population_string += ("\n\twriteFile(\"" + population_parameters["pop_name"] + "_fixed_mutation_counts.txt\"," +
                "asString(sim.getValue(\"fixations_counted_p1\")));" +
                "\n\twriteFile(\"" + population_parameters["pop_name"] + "_fixed_mutations.txt\"," +
                "\n\tpaste(sim.getValue(\"fixations_p1\"), sep = \"\"));")

        end_population_string += "\n\tsim.outputFixedMutations();"


        #If this is the last clade from a certain parent, write script to destroy that parent's temporary files
        if(population_parameters["last_child_clade"]):
            end_population_string += ("\n\n\tsystem(\"rm " + population_parameters["parent_pop_name"] + ".txt\");" +
                                      "\n\tsystem(\"rm " + population_parameters["parent_pop_name"] + ".fasta\");")


        if(population_parameters["terminal_clade"]):
            end_population_string += "\n}"

        self.output_file.write(end_population_string)
