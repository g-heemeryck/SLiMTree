

#Program to write SLiM Scripts from set functions
#Required packages:
#random
#csv
#numpy
#os

import random, csv, os
import numpy as np

#Base class from which other classes inherit
class writeSLiM:

    #Initialize required parameters
    def __init__(self, start_para_dict, partition_information, repeated_commands_booleans):

        #Set up variables that remain constant for every part of the simulation
        self.general_output_filename = start_para_dict["output_file"]
        self.genome_length = start_para_dict["genome_length"]
        self.sample_size = start_para_dict["sample_size"]
        self.fasta_filename = start_para_dict["fasta_filename"]
        self.partition = partition_information[0]
        self.partition_time = partition_information[1]
        self.repeated_commands_booleans = repeated_commands_booleans
        
        #Set up the fitness profile and starting distribution of amino acids
        self.fitness_profile_nums = start_para_dict["fitness_profile_nums"]
        self.fitness_profiles = start_para_dict["fitness_profiles"]
        self.starting_allele_dist = start_para_dict["stationary_distributions"]
        self.amino_acids = start_para_dict["amino_acids"]
        
        #Set up type of model?
        self.model_type = start_para_dict["wf_model"]

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



        #Write the initialize function for the SLiM script
    def write_initialize(self, population_parameters):
        
        initialize_string = ("initialize() {")

        if (self.model_type == False): #***
            initialize_string += "\n\tinitializeSLiMModelType(\"nonWF\");" #****

        initialize_string += ("\n\tsetSeed(" + str(random.randint(0,1000000000)) + ");" + "\n\tinitializeSLiMOptions(nucleotideBased=T);")

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
                           "\n\tdefineConstant(\"seq_length\", " + str(len(self.fitness_profile_nums)+2) + ");")

        #Write code to start a fixed state from the starting nucleotide sequence
        set_up_fitness += "sim.setValue(\"fixations_p1\", strsplit(sim.chromosome.ancestralNucleotides(),sep = \"\"));"

        #At the start of the sim there are no fixations counted
        set_up_fitness += "\n\tsim.setValue(\"fixations_counted_p1\", 0);"
        set_up_fitness += "\n}\n\n\n"

                                                                                                  
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


    #Write the reproduction callback for non-Wright-Fisher models
    def write_reproduction(self):
        #Basic reproduction callback for now; more functionality could be added later if necessary.
        reproduction_string = ("reproduction() { " +
                            "\n\tsubpop.addCrossed(individual, subpop.sampleIndividuals(1));" +
                            "\n }\n\n\n")

        self.output_file.write(reproduction_string)

    #Write code to count substitutions, make a backup and count generations
    def write_repeated_commands(self, start_dist, end_dist, pop_name, count_subs = True, output_gens = True, backup = True):
        repeated_commands_string = str(start_dist) +":" + str(end_dist) + "late () {"
        
        #Write a command to count the substitutions (identity by state)
        if (count_subs):
            repeated_commands_string += ("\n\tif(length(sim.mutations)!= 0){"
                        "\n\t\tancestral_genome = sim.getValue(\"fixations_" + pop_name + "\");" +
                        "\n\t\tcompare_genome = strsplit(" + pop_name + ".genomes[0].nucleotides(), sep = \'\');"+
                        "\n\t\tfixed_nucs = rep(T, length(compare_genome));" +
                        "\n\n\t\tfor (genome in (" + pop_name + ".genomes)){" +
                        "\n\t\t\tsame_nucs = (compare_genome == strsplit(genome.nucleotides(), sep = \'\'));" +
                        "\n\t\t\tfixed_nucs = (fixed_nucs & same_nucs);\n\t\t}" +
                        "\n\n\t\tdifferent_muts = (ancestral_genome != compare_genome);" +
                        "\n\t\tnew_fixations = different_muts & fixed_nucs;" +
                        "\n\t\tsim.setValue(\"fixations_counted_" + pop_name + 
                        "\", sim.getValue(\"fixations_counted_" + pop_name+ "\") + sum(new_fixations));" +
                        "\n\n\t\tancestral_genome[new_fixations] = compare_genome[new_fixations];" +
                        "\n\t\tsim.setValue(\"fixations_" + pop_name + "\", ancestral_genome);\n\t};")
          
        #Write a command to output when every 100th generation has passed
        if(output_gens):
            repeated_commands_string += "\n\n\tif (sim.generation%100 == 0) {\n\t\tcatn(sim.generation);\n\t};"
            

        #Write a command to write a backup of all individuals after every 100 generations
        if (backup):
             repeated_commands_string += ("\n\n\tif (sim.generation%100 == 0) {" +
                        "\n\t\twriteFile(\"" + os.getcwd()+ "/backupFiles/" + pop_name + ".fasta\"," + 
                        "(\">parent_ancestral_to_load\\n\" + sim.chromosome.ancestralNucleotides()));" +
                        "\n\t\tsim.outputFull(\"" + os.getcwd()+ "/backupFiles/" + pop_name + ".txt\");\n\t};")
        
        repeated_commands_string += "\n}\n\n\n"
        
        self.output_file.write(repeated_commands_string)
        

    
    #Write code to add first population, subpopulation or completely remove population and replace with another
    def write_subpop(self, population_parameters):
            
                                        
        if(population_parameters["parent_pop_name"] == None):
                self.set_up_sim(population_parameters)
        else:
            #Not the starting population, break off from existing population
            define_population_string = (str(int(population_parameters["dist_from_start"])+1) + " { \n" +
                    "\tsim.addSubpopSplit(\""+ population_parameters["pop_name"] + "\"," +
                    str(population_parameters["population_size"]) + ", " + population_parameters["parent_pop_name"]+ ");"+
                    "\n\n\tsim.setValue(\"fixations_" + population_parameters["pop_name"] + "\", sim.getValue(\"fixations_"+
                    population_parameters["parent_pop_name"] +"\"));" +
                    "\n\tsim.setValue(\"fixations_counted_"+ population_parameters["pop_name"]+"\", 0);")

            if(population_parameters["last_child_clade"] == True):
                define_population_string += "\n\t" + population_parameters["parent_pop_name"]+".setSubpopulationSize(0);"    

            define_population_string += "\n}\n\n\n" 
                                            
            self.output_file.write(define_population_string)
                
                
        #Write the commands that are run for every simulation and the starting population
        self.write_repeated_commands(int(population_parameters["dist_from_start"])+1, 
                        int(population_parameters["end_dist"]), population_parameters["pop_name"],
                        self.repeated_commands_booleans[0], self.repeated_commands_booleans[1],
                        self.repeated_commands_booleans[2])
                        
        #Write the end of each population
        self.write_end_pop(population_parameters)
        
    #Write code to add first population, subpopulation or completely remove population and replace with another with non-Wright-Fisher models
    def write_subpop_nonwf(self, population_parameters):
        if(population_parameters["parent_pop_name"] == None):
                self.set_up_sim(population_parameters)
        else:
            #Not the starting population, break off from existing population
            define_population_string = (str(int(population_parameters["dist_from_start"])) + " late() { \n" +
                                    "\tsim.addSubpop(\"" + population_parameters["pop_name"] + "\", 0);")
            #If this is the last population broken off, take the last half of the parent population
            if (population_parameters["last_child_clade"] == True):
                define_population_string += str("\n\t" + population_parameters["pop_name"] + ".takeMigrants(" + population_parameters["parent_pop_name"] + ".individuals);")
            else:
                #Take half of the parent population
                define_population_string += str("\n\tmigrants = sample(" + population_parameters["parent_pop_name"] + ".individuals, integerDiv("
                                    + population_parameters["parent_pop_name"] + ".individualCount, 2));\n\t"
                                    + population_parameters["pop_name"] + ".takeMigrants(migrants);")

            define_population_string += str("\n\n\tsim.setValue(\"fixations_" + population_parameters["pop_name"] + "\", sim.getValue(\"fixations_"+
                                    population_parameters["parent_pop_name"] +"\"));" +
                                    "\n\tsim.setValue(\"fixations_counted_"+ population_parameters["pop_name"]+"\", 0);")


            define_population_string += "\n}\n\n\n"

            self.output_file.write(define_population_string)

        #Write the early commands - this may need tweaking w/ the fitness algorithm

        early_event = str(int(population_parameters["dist_from_start"]) + 1) + ":" + str(int(population_parameters["end_dist"])) + " early(){\n\t" + population_parameters["pop_name"] + ".fitnessScaling = 500/ " + population_parameters["pop_name"] + ".individualCount;" + "\n}\n\n\n"

        self.output_file.write(early_event)


        #Write the commands that are run for every simulation and the starting population
        self.write_repeated_commands(int(population_parameters["dist_from_start"])+1,
                        int(population_parameters["end_dist"]), population_parameters["pop_name"],
                        self.repeated_commands_booleans[0], self.repeated_commands_booleans[1],
                        self.repeated_commands_booleans[2])

        #Write the end of each population
        self.write_end_pop(population_parameters)
    
    
    #Set up the simulation by initializing everything
    def set_up_sim(self, population_parameters):
        self.output_file = open(self.general_output_filename + "_" + population_parameters["pop_name"] + ".slim" , "w")

        #Set up the initialize and fitness functions for the new script
        self.write_initialize(population_parameters)
        self.write_fitness()
        
        #Write reproduction callback if this is a non-WF model
        if (self.model_type == False):
            self.write_reproduction()
        
        #Make the population and set up fitness effects
        pop_string = ("1 early() {" +
                    "\n\tsetup_fitness();" +
                    "\n\twriteFile(\"" + self.fasta_filename + "_aa.fasta\", \"\", append = F);" +
                    "\n\twriteFile(\"" + self.fasta_filename + "_nuc.fasta\", \"\", append = F);" +
                    "\n\tsim.addSubpop(\"p1\", " + str(population_parameters["population_size"]) + ");")
        
        #Write code to start a fixed state from the starting nucleotide sequence
        pop_string += "\n\tsim.setValue(\"fixations_p1\", strsplit(sim.chromosome.ancestralNucleotides(),sep = \"\"));"

        #At the start of the sim there are no fixations counted
        pop_string += "\n\tsim.setValue(\"fixations_counted_p1\", 0);"
        pop_string += "\n}\n\n\n"
                                            
        self.output_file.write(pop_string)
        
        
    #Write the end of a population to save the number of substitutions and output sequence data
    def write_end_pop (self, population_parameters):
        end_population_string = str(int(population_parameters["end_dist"])) + " late() {"

        #If terminal clade output data 
        if(population_parameters["terminal_clade"]):
            end_population_string += self.write_terminal_output(population_parameters, pop = population_parameters["pop_name"])
        
        #Write file with the substitution counts
        if(self.repeated_commands_booleans[0]):
            end_population_string += ("\n\twriteFile(\"" + os.getcwd()+ "/" + population_parameters["pop_name"] + "_fixed_mutation_counts.txt\"," +
                "asString(sim.getValue(\"fixations_counted_" + population_parameters["pop_name"] + "\")));" +
                "\n\twriteFile(\"" + os.getcwd()+ "/" + population_parameters["pop_name"] + "_fixed_mutations.txt\"," +
                " paste(sim.getValue(\"fixations_" + population_parameters["pop_name"] + "\"), sep = \"\"));")
                
        end_population_string += "\n}\n\n\n"        
        
        self.output_file.write(end_population_string)
 



    #Write code to write the output for terminal populations after they have reached their population
    def write_terminal_output(self, population_parameters, pop = "p1"):

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
            
            
            terminal_output_string += ("\n\tg = "+ pop +".genomes[" + str(sample) + "];" +
                                       "\n\twriteFile(\"" + nuc_filename + "\", " + fasta_string_nuc +", append = T);" +
                                       "\n\twriteFile(\"" + aa_filename + "\", " + fasta_string_aa +", append = T);\n\n")
            count += 1


        return terminal_output_string	

