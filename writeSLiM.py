

#Program to write SLiM Scripts from set functions
#Required packages:
#random
#csv
#numpy
#os

import random, csv, os
import numpy as np

class writeSLiM:

        #Initialize required parameters
        def __init__(self, start_para_dict, partition_information):

                #Set up variables that remain constant for every part of the simulation
                self.general_output_filename = start_para_dict["output_file"]
                self.genome_length = start_para_dict["genome_length"]
                self.sample_size = start_para_dict["sample_size"]
                self.fasta_filename = start_para_dict["fasta_filename"]
                self.partition = partition_information[0]
                self.partition_time = partition_information[1]
                
                #Set up the fitness profile and starting distribution of amino acids
                self.fitness_profile_nums = start_para_dict["fitness_profile_nums"]
                self.fitness_profiles = start_para_dict["fitness_profiles"]
                self.starting_allele_dist = start_para_dict["stationary_distributions"]
                self.amino_acids = start_para_dict["amino_acids"]


                #Set up conversion from amino acid to codon
                with open(os.path.dirname(os.path.realpath(__file__)) + '/fitnessDataFiles/slim_codon_nums.csv', newline='') as slim_codon_nums:
                        reader = csv.reader(slim_codon_nums)
                        slim_codons = list(reader)[1:]

                self.slim_codon_dict = {}

                for codons in slim_codons :
                        amino_acid = codons[2]
                        slim_codon_number = int(codons[0])
                        
                        if(amino_acid in self.slim_codon_dict.keys()):
                                self.slim_codon_dict[amino_acid].append(slim_codon_number)
                        else:
                                self.slim_codon_dict[amino_acid] = [slim_codon_number]





        #Write code to add a new population to the simulation by writing a new script for that population
        def write_subpop(self, population_parameters):
                
                #Create a new script and batch file for the population
                batch_file = open(self.general_output_filename + "_" + population_parameters["pop_name"] + ".sh", "w")
                batch_file.write("#!/bin/sh\n\n#SBATCH -J SLiM_Simulation_" + population_parameters["pop_name"] + "\n#SBATCH -t " + self.partition_time + 
                        "\n#SBATCH -p "  + self.partition + "\n#SBATCH -o " + population_parameters["pop_name"] + ".out\n#SBATCH -e " + population_parameters["pop_name"] 
                        +".err\n#SBATCH -n 1" + "\n\nslim " + self.general_output_filename + "_" + population_parameters["pop_name"]+".slim")
                batch_file.close()

                self.output_file = open(self.general_output_filename + "_" + population_parameters["pop_name"] + ".slim" , "w")

                #Set up the initialize and fitness functions for the new script
                self.write_initialize(population_parameters)
                self.write_fitness()

                #Write the commands to count the substitutions, track the generation and save a backup file
                simulation_distance_string = (str(int(population_parameters["dist_from_start"])+2) +":" + str(int(population_parameters["end_dist"])) +
                        "late (){\n\tunique_mutations_bools = sapply(sim.mutations, \"sum(applyValue.position == sim.mutations.position);\") == 1;" +
                        "\n\tancestral_genome = sim.getValue(\"fixations\");" +
                        "\n\tcompare_genome = strsplit(p1.genomes[0].nucleotides(), sep = \'\');"+
                        "\n\tfixed_nucs = rep(T, length(compare_genome));" +
                        "\n\n\tfor (genome in (p1.genomes)){" +
                        "\n\t\tsame_nucs = (compare_genome == strsplit(genome.nucleotides(), sep = \'\'));" +
                        "\n\t\tfixed_nucs = (fixed_nucs & same_nucs);\n}" +
                        "\n\n\tdifferent_muts = (ancestral_genome != compare_genome);" +
                        "\n\tnew_fixations = different_muts & fixed_nucs;" +
                        "\n\tsim.setValue(\"fixations_counted\", sim.getValue(\"fixations_counted\") + sum(new_fixations));" +
                        "\n\n\tancestral_genome[new_fixations] = compare_genome[new_fixations];" +
                        "\n\tsim.setValue(\"fixations\", ancestral_genome);"+
                        "\n\n\tif (sim.generation%100 == 0) {\n\t\tcatn(sim.generation);"+
                        "\n\t\twriteFile(\"backupFiles/" + population_parameters["pop_name"] +
                        ".fasta\", (\">parent_ancestral_to_load\\n\" + sim.chromosome.ancestralNucleotides()));" +
                        "\n\t\tsim.outputFull(\"backupFiles/" + population_parameters["pop_name"] + ".txt\");\n\t};\n}\n\n\n")
                self.output_file.write(simulation_distance_string)
                self.write_start_pop(population_parameters)


                #Finish writing the script
                self.write_end_sim(population_parameters)
                self.output_file.close() 

                




        #Write the initialize function for the SLiM script
        def write_initialize(self, population_parameters):
                
                initialize_string = ("initialize() {" +
                                     "\n\tsetSeed(" + str(random.randint(0,1000000000)) + ");" +
                                     "\n\tinitializeSLiMOptions(nucleotideBased=T);")


                #Starting population does not inherit parent sequence, other populations do
                if(population_parameters["parent_pop_name"] == None):                
                        aa_codon_sequence = str(self.create_codon_seq())
                        aa_codon_sequence_str = "c(" + aa_codon_sequence[1:len(aa_codon_sequence)-1] + ")"
                        initialize_string += ("\n\tcodons = " + aa_codon_sequence_str + ";" +
                                        "\n\tinitializeAncestralNucleotides(codonsToNucleotides(codons, format=\"char\"));" +
                                        "\n\tdefineConstant(\"L\"," + str(int(self.genome_length*3)) + ");")
                else:
                        initialize_string += ("\n\tdefineConstant(\"L\", initializeAncestralNucleotides(\"" +
                                              population_parameters["parent_pop_name"] + ".fasta\"));")

                initialize_string += ("\n\tmm = mmJukesCantor(" + str(population_parameters ["mutation_rate"]/3) + ");" +
                                "\n\tinitializeMutationTypeNuc(\"m1\", 0.5, \"f\", 0.0);" +
                                "\n\tm1.convertToSubstitution = F;" + 
                                "\n\tinitializeGenomicElementType(\"g1\", m1, 1.0, mm);" +
                                "\n\tinitializeGenomicElement(g1, 0, L-1);" +
                                "\n\tinitializeRecombinationRate("+ str(population_parameters ["recombination_rate"])+");"+
                                "\n}\n\n\n")
                
                self.output_file.write(initialize_string)

        



        #Create an initial codon sequence to put into SLiM based on the fitness profile
        def create_codon_seq(self):
                
                #Methionine - start codon
                start_codon = [14]

                #Stop codons
                stop_codons = [48, 50, 56]
                stop_codon = [random.choice(stop_codons)]

                #Middle codons - chosen according to distribution of alleles
                middle_amino_acids = []

                for dist_num in self.fitness_profile_nums:
                        weights = self.starting_allele_dist[:,dist_num].tolist()
                        middle_amino_acids += random.choices(self.amino_acids, weights = weights, k = 1)

                middle_codons = list(map(self.convert_amino_acid, middle_amino_acids))


                aa_sequence_in_codons = start_codon + middle_codons + stop_codon
                return (aa_sequence_in_codons)




        #Convert an amino acid to a codon in SLiM by choosing a random available codon
        def convert_amino_acid(self, amino_acid):
                codon_list = self.slim_codon_dict[amino_acid]
                selected_codon = codon_list[random.randint(0,len(codon_list)-1)]

                return selected_codon



        #Write the fitness callback for the SLiM script according to the distribution of fitness effects
        def write_fitness(self):
                #Set up a dictionary in SLiM which takes in amino acids as keys and returns vector of fitnesses
                set_up_fitness = "function (void) setup_fitness(void){\n"

                for key_value in self.fitness_profiles:
                        aa_fitnesses = str(self.fitness_profiles[key_value])
                        aa_fitnesses = "c(" + aa_fitnesses[1:len(aa_fitnesses)-1] + ")"
                        set_up_fitness += "\tsim.setValue(\"" + key_value + "\", " + aa_fitnesses + ");\n"

                fitness_vector = str(self.fitness_profile_nums)
                fitness_vector = "c(" + fitness_vector[1:len(fitness_vector)-1] + ")"
                set_up_fitness += ("\n\tdefineConstant(\"fitness_profiles\", " + fitness_vector +");" +
                                   "\n\tdefineConstant(\"seq_length\", " + str(len(self.fitness_profile_nums)+2) + ");" +
                                   "\n}\n\n\n")
                                                                                                          
                self.output_file.write(set_up_fitness)

                #Sets up a function to return the fitness of each individual amino acid which scales stop codons by position
                get_aa_fitness = ("function (float) get_aa_fitness (string aa_seq, integer position){" +
                                "\n\tif (aa_seq [position] == \"X\"){"+
                                "\n\t\treturn position*(1/seq_length);" +
                                "\n\t} else {" +
                                "\n\t\treturn sim.getValue(aa_seq[position])[fitness_profiles[position]];" +
                                "\n\t}\n}\n\n\n")
                self.output_file.write(get_aa_fitness)


                
                #Defining a function in SLiM which returns the fitness of the amino acid sequence
                fitness_function_string = ("function (float) get_fitness (string aa_seq_string){" +
                        "\n\taa_seq = strsplit(aa_seq_string, sep=\"\");" +
                        "\n\n\tif(aa_seq[0] != \"M\" | aa_seq[seq_length-1] != \"X\"){" +
                        "\n\t\treturn 0.1;\n\t}" +
                        "\n\n\taa_seq = aa_seq[1:(seq_length-2)];" +
                        "\n\tfitnesses = sapply(seq(0,seq_length-3), " +
                                           "\"sim.getValue(aa_seq[applyValue])[fitness_profiles[applyValue]];\");" +
                        "\n\n\treturn product(fitnesses);\n}\n\n\n")
                        
                self.output_file.write(fitness_function_string)


                
               #Now write out the fitness callback based on the fitness distribution
                fitness_callback_string = ("fitness(NULL) {" +
                                "\n\tfor (g in individual.genomes){" +
                                "\n\t\taa_seq = codonsToAminoAcids(nucleotidesToCodons(g.nucleotides()));"+
                                "\n\t\treturn get_fitness(aa_seq);" +
                                "\n\t}\n}\n\n\n")
                
                self.output_file.write(fitness_callback_string)




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
                pop_string += "sim.setValue(\"fixations_counted\", 0);"


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









		
		
	
