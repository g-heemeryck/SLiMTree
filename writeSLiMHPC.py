#Program to write SLiMTree code for a high performance computing cluster


import random, csv, os
from writeSLiM import writeSLiM



class writeSLiMHPC(writeSLiM):
    #Write code to add a new population to the simulation by writing a new script for that population
    def write_subpop(self, population_parameters, ignore):
        
        #Create a new script and batch file for the population
        batch_file = open(self.general_output_filename + "_" + population_parameters["pop_name"] + ".sh", "w")
        batch_file.write("#!/bin/sh\n\n#SBATCH -J SLiM_Simulation_" + population_parameters["pop_name"] + "\n#SBATCH -t " + self.partition_time + 
                "\n#SBATCH -p "  + self.partition + "\n#SBATCH -o " + population_parameters["pop_name"] + ".out\n#SBATCH -e " + population_parameters["pop_name"] 
                +".err\n#SBATCH -n 1" + "\n\nslim " + self.general_output_filename + "_" + population_parameters["pop_name"]+".slim")
        batch_file.close()

        self.output_file = open(self.general_output_filename + "_" + population_parameters["pop_name"] + ".slim" , "w")

        #Set up the initialize and fitness functions for the new script
        super().write_initialize(population_parameters)
        super().write_fitness()

        #Write the commands that are run for every simulation and the starting population
        super().write_repeated_commands(int(population_parameters["dist_from_start"])+2, 
                            int(population_parameters["end_dist"]), population_parameters["pop_name"],
                            self.repeated_commands_booleans[0], self.repeated_commands_booleans[1], 
                            self.repeated_commands_booleans[2])
        self.write_start_pop(population_parameters)


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
            pop_string += "sim.setValue(\"fixations\", strsplit(sim.chromosome.ancestralNucleotides(),sep = \"\"));"
        else:
            pop_string += ("\n\tsim.readFromPopulationFile(\"" + population_parameters["parent_pop_name"]  + ".txt\");" +
                            "\n\tp1.setSubpopulationSize(" + str(population_parameters["population_size"]) + ");")

            #Load population into the end of the parent population's script to start this script when parent's finishes
            parent_output_file = open(self.general_output_filename + "_" + population_parameters["parent_pop_name"] + ".slim" , "a")
            parent_output_file.write("\n\tsystem(\"sbatch \\\"" + self.general_output_filename + "_" + population_parameters["pop_name"] + ".sh\\\"\");")

            if(population_parameters["last_child_clade"]):
                parent_output_file.write("\n}")

            #Write code to import in the prevouisly fixed state
            pop_string += ("sim.setValue(\"fixations\", strsplit(readFile(\""+ population_parameters["parent_pop_name"] +
                           "_fixed_mutations.txt\"), sep = \"\"));")

        #At the start of the sim there are no fixations counted
        pop_string += "\n\tsim.setValue(\"fixations_counted\", 0);"
        pop_string += "\n}\n\n\n"
                                            
        self.output_file.write(pop_string)
        


    #Write the closing statements to end the simulation and either allow for starting of subsequent simulations or output data (depending on terminal status)
    def write_end_sim(self, population_parameters):
        
        end_population_string = str(int(population_parameters["end_dist"]) + 1) + " late() {"

        #If terminal clade output data otherwise create data to be loaded into scripts of the clades children
        if(population_parameters["terminal_clade"]):
            end_population_string += self.write_terminal_output(population_parameters)
        else:
            end_population_string += ("\n\twriteFile(\"" + population_parameters["pop_name"] +
                                      ".fasta\", (\">parent_ancestral_to_load\\n\" + sim.chromosome.ancestralNucleotides()));" +
                                      "\n\tsim.outputFull(\"" + population_parameters["pop_name"] + ".txt\");" )

        #Scripting to end the simulation and write the fixed mutations
        end_population_string += ("\n\twriteFile(\"" + population_parameters["pop_name"] + "_fixed_mutation_counts.txt\"," +
                "asString(sim.getValue(\"fixations_counted\")));" +
                "\n\twriteFile(\"" + population_parameters["pop_name"] + "_fixed_mutations.txt\"," +
                " paste(sim.getValue(\"fixations\"), sep = \"\"));")

        end_population_string += "\n\tsim.outputFixedMutations();"


        #If this is the last clade from a certain parent, write script to destroy that parent's temporary files
        if(population_parameters["last_child_clade"]):
            end_population_string += ("\n\n\tsystem(\"rm " + population_parameters["parent_pop_name"] + ".txt\");" +
                                      "\n\tsystem(\"rm " + population_parameters["parent_pop_name"] + ".fasta\");")

        
        if(population_parameters["terminal_clade"]):
            end_population_string += "\n}"
        
        self.output_file.write(end_population_string)



    #Write code to write the output for terminal populations after they have reached their population
    def write_terminal_output(self, population_parameters):

        #Set up the names of the 3 fasta files to be output to
        nuc_filename = self.fasta_filename + "_nuc.fasta"
        aa_filename =  self.fasta_filename + "_aa.fasta"
        ancestral_filename = self.fasta_filename + "_fixed.fasta"


        #Set up sampling of the population
        pop_name = population_parameters["pop_name"]              
        pop_size = population_parameters["population_size"]

        if(self.sample_size == "all"):
            pop_samples = list(range(pop_size))
        elif(int(self.sample_size) < pop_size):
            pop_samples = random.sample(list(range(pop_size)), int(self.sample_size))
        else:
            pop_samples = list(range(pop_size))



        #Iterate through each random sample to write script to output samples of amino acids and nucleotides to fasta files
        terminal_output_string = ""
        count = 0
        for sample in pop_samples:
            fasta_string_nuc = "\">" + pop_name + "_" + str(count) + ": \\n\" + g.nucleotides()"
            fasta_string_aa = "\">" + pop_name + "_" + str(count) + ": \\n\" + codonsToAminoAcids(nucleotidesToCodons(g.nucleotides()))"
            
            
            terminal_output_string += ("\n\tg = p1.genomes[" + str(sample) + "];" +
                                       "\n\twriteFile(\"" + nuc_filename + "\", " + fasta_string_nuc +", append = T);" +
                                       "\n\twriteFile(\"" + aa_filename + "\", " + fasta_string_aa +", append = T);\n\n")
            count += 1


        return terminal_output_string	