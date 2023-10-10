def parseConfigFile(filename):
    config_dictionary = {}
    # create a dict to store parsed strings
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            # skip over uneeded lines
            if line.startswith(('#', ';')) or not line:
                continue
            if '=' in line:
                key, value = [item.strip() for item in line.split('=')]
                config_dictionary[key] = value
                # create key value pairs
    return config_dictionary

def readFASTA(filename):
    sequences = []
    with open(filename, 'r') as file:
        seq_name, description, sequence = None, None, ''
        # initialize important strings
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if seq_name and description and sequence:  # Save the previous sequence
                    sequences.append((seq_name, description, sequence))
                    # if line starts with > append to list
                seq_name, description = line[1:].split(' ', 1)  # Remove ">" and split at the first space
                sequence = ''  # Initialize the sequence
            else:
                sequence += line
        if seq_name and description and sequence:  # Save the last sequence
            sequences.append((seq_name, description, sequence))
    return sequences


def printScreen1(config, sequences):
    print("Welcome Sequence Viewer!")
    print(f"Programmer: {config['Programmer']}")
    print(f"Email: {config['Email']}\n")
    print(f"There is a total of {len(sequences)} sequences detected in the file: {config['SeqFileName']}\n")
    # print out welcome mesage for screen1


def getNucleotidesPerLine(config_dict):
    return int(config_dict['NucleotidesPerLine[50|100]'])
    # check how many nucleotides to display per line

def printInFASTA(seq_name, sequence, description, NucleotidesPerLine):
    print(f">{seq_name} {description}")
    for i in range(0, len(sequence), NucleotidesPerLine):
        print(sequence[i:i+NucleotidesPerLine])
        # print nucleotides according to how many per line in config file
    print()

def printWithRuler(seq_name, Sequence, description, Spacer, NucleotidesPerLine):
    seq_length = len(Sequence)
    num_lines = (seq_length + NucleotidesPerLine - 1) // NucleotidesPerLine
    seq_name = seq_name
    description = description

    # print the name and description
    print(f">{seq_name} {description}\n")

    # print first ruler with numbers every ten nucleotides
    ruler1 = " " * 6
    for i in range(1, 11):
        ruler1 += f"{i:10}{Spacer}"  # add the spacer
    print(ruler1)

    # print second ruler with index numbers
    ruler2 = ""
    for i in range(1, 101):
        ruler2 += f"{i % 10}"
        if i % 10 == 0:
            ruler2 += Spacer  # spacer every ten numbers
    print("Line  " + ruler2)

    # print the sequence lines
    for line_number in range(num_lines):
        start = line_number * NucleotidesPerLine
        end = min(start + NucleotidesPerLine, seq_length)
        line = Sequence[start:end]

        # spacer added every 10
        line_with_spacer = []
        for i in range(0, len(line), 10):
            line_with_spacer.append(line[i:i + 10])
        line_number_str = str(line_number + 1).rjust(4)

        # print the line with spacer
        print(f" {line_number_str} {Spacer.join(line_with_spacer)}")

    print()


def nucleotideCounter(sequence):
    return (len(sequence), sequence.count('A'), sequence.count('T'), sequence.count('G'), sequence.count('C'), sequence.count('N'))
    # count how many of each nucleotide and return to be printed

def gcContent(Sequence):
    gc = Sequence.count('G') + Sequence.count('C')
    total = len(Sequence)
    return round((gc/total) * 100, 2)
    # calculate gc content

def diNucleotideProfile(sequence):
    dinucleotides = ['AA', 'AT', 'AG', 'AC', 'TA', 'TT', 'TG', 'TC', 'GA', 'GT', 'GG', 'GC', 'CA', 'CT', 'CG', 'CC']
    profile = {key: sequence.count(key) for key in dinucleotides}
    # create key for 16 dinucleotides

    formatted_profile = " ".join([f"{key}={profile[key]}" for key in dinucleotides])
    # return in correct format

    return formatted_profile



def CpGIsland(sequence):
    islands = {}
    cpg = 'CG' # looking for CG in sequence
    start = 0
    island_count = 0 # counting number of CG
    formatted_islands = ""

    while start < len(sequence) - 6:
        ind = sequence.find(cpg, start) # where CG first is seen
        if ind == -1:
            break
        end = ind
        while sequence[end:end + 2] == cpg:
            end += 2
        if end - ind >= 6:
            island_count += 1
            island_info = f"{ind}-{end - 1}_{end - ind}"
            islands[island_count] = island_info
            formatted_islands += f"{island_count}={island_info} "
            start = end
        else:
            start = ind + 2

    return formatted_islands

def processInquiry(config, sequences):
    selected_seqs = [int(ind) for ind in config['SelectedSeqs'].split(',')]
    start_position = int(config['SeqFragmentStartPosition'])
    end_position = int(config['SeqFragmentEndPostion'])
    # get the fragment based on the positions in config file
    total_seqs = len(sequences)
    selected_seqs = [ind for ind in selected_seqs if 1 <= ind <= total_seqs]

    if not selected_seqs:
        print("No valid sequences selected for inquiry.")
        return

    print(f"Among the {total_seqs} sequences detected in the file: {config['SeqFileName']}")
    print(f"You have selected {selected_seqs} for the inquiry mode.")
    print(f"The start and end positions for sequence fragments: {start_position}-{end_position}\n")

    # print out the message for user

    for ind in selected_seqs:
        seq_name, description, sequence = sequences[ind - 1]

        # get the fragment
        selected_fragment = sequence[start_position - 1:end_position]

        print(f">{seq_name} {description}\n")

        # print the fragment with ruler
        ruler = f"<{start_position}{'-' * (end_position - start_position + 1)}{end_position}>"
        print(ruler)
        print("|" * (end_position - start_position + 2))
        print(selected_fragment)

        # print nucleotide counts
        count = nucleotideCounter(selected_fragment)
        print(f"Nucleotide Counts: Seq Length={count[0]} A={count[1]} T={count[2]} G={count[3]} C={count[4]} N={count[5]}")

        # print gc content
        print(f"GC Content={gcContent(selected_fragment)}%")

        # print dinucleotide profile
        print(f"Dinucleotide profile: {diNucleotideProfile(selected_fragment)}")

        # print cpG islands
        print(f"CpG Islands: {CpGIsland(selected_fragment)}")

        print()



if __name__ == '__main__':
    import sys

    config_filename = sys.argv[1]
    config = parseConfigFile(config_filename)
    sequences = readFASTA(config["SeqFileName"])
    NucleotidesPerLine = getNucleotidesPerLine(config)

    # print screen 1
    printScreen1(config, sequences)

    for seq_name, description, sequence in sequences:
        # loop through all sequences
        if config["ViewSequenceInFastaFormat[Y|N]"] == 'N': # print based on ruler needed
            if config["DoYouNeedSpaceSeperator[Y|N]"] == 'Y': # print based on spacer needed
                spacer = ' '
            else:
                spacer = ''
            printWithRuler(seq_name, sequence, description, spacer, NucleotidesPerLine)
        else:
            printInFASTA(seq_name, sequence, description, NucleotidesPerLine)

        # print nucleotide count
        if config["nucleotideCounter[Y|N]"] == 'Y':
            count = nucleotideCounter(sequence)
            print(
                f"Nucleotide Counts: Seq Length={count[0]} A={count[1]} T={count[2]} G={count[3]} C={count[4]} N={count[5]}")

        # print gc content
        if config["gcContent[Y|N]"] == 'Y':
            print(f"GC content={gcContent(sequence)}%")

        # print dinucleotide profile
        if config["dinucleotideProfil[Y|N]"] == 'Y':
            print(f"Dinucleotide profile: {diNucleotideProfile(sequence)}")

        #  print cpg islants
        if config["CpGIsland[Y|N]"] == 'Y':
            print(f"CpG Islands: {CpGIsland(sequence)}")

        print()

    processInquiry(config, sequences)
    # run the inquiry mode