import sys
from Bio import SeqIO
import re
import editdistance
import random
import argparse
from itertools import combinations

def main():
    parser = argparse.ArgumentParser(prog="disco",add_help=False)
    parser.add_argument("-i", "--input", type=argparse.FileType("r"),dest="input",
                      help="alignment file in fasta form")
    parser.add_argument("--editdistance",type=int, dest="edit_value",
                     help="Edit distance value as interger")
    parser.add_argument("-o","--output",type=argparse.FileType("w"),dest="output",
                        help="output file name")
    parser.add_argument("--seed",type=int,default=10,
                        help="seed number for reproducibility, default is 10")


    subparsers = parser.add_subparsers(help="sub-command help")

    #Create subcommand
    parser_create = subparsers.add_parser("create", parents=[parser],help='Module to create highly diverse community at specified edit distance')
    parser_create.add_argument("--community",type=argparse.FileType("r"),dest="starter_community",
                        help="Starting community with each identifier on its own line")
    parser_create.add_argument("--trimPrimers",type=int,dest="trimPrimers",
                        help="length of primers to trim")
    parser_create.add_argument("--metadata",type=argparse.FileType("r"),
                        help="information to combine with the community output, file must contain information in the first column and the identifiers in the last in tab deliminated form, with a header")
    parser_create.add_argument("--fasta",type=argparse.FileType("w"),dest="output_fasta",
                        help="Output final community fasta file")
    parser_create.set_defaults(func=create)

    #Correct subcommand
    parser_correct=subparsers.add_parser("correct", parents=[parser],help="Module to correct sequencing errors based on the community and edit distance")
    parser_correct.add_argument("--EDC",type=argparse.FileType("r"),dest="input_community",
                        help="Input community with each identifier on its own line")
    parser_correct.set_defaults(func=correct)

    #Subsample subcommand
    parser_subsample=subparsers.add_parser("subsample", parents=[parser],help="Module to subsample highly diverse community")
    parser_subsample.add_argument("--strains",type=int, dest="strains",
                    help="number of strains desired in final community")
    parser_subsample.add_argument("--taxonomy",type=argparse.FileType("r"),dest="taxonomy",
                        help="Tab seperated file with strain ids in first column and taxonomic rank in second")
    parser_subsample.add_argument("--proportion",type=argparse.FileType("r"),dest="proportion",
                        help="File of the relative proportions of each taxonomic rank desired in final community")
    parser_subsample.set_defaults(func=subsample)

    # Parse args
    args = parser.parse_args()

def create(args):
    random.seed(args.seed)
    if (args.input):
        if (args.trimPrimers):
            print("Creating sequence dictionary")
            sequence_dict=TrimPrimers(args.input,args.trimPrimers)
            print("Creating edit distance dictionary")
            dict_ed=editDistanceDictionary(sequence_dict)
        else:
            print("Creating sequence dictionary")
            sequence_dict=sequenceDictionary(args.input)
            print("Creating edit distance dictionary")
            dict_ed=editDistanceDictionary(sequence_dict)
    else:
        sys.exit("ERROR:No input file")
    if (args.edit_value):
        if (args.starter_community):
            starter_community=startcommunity(args.starter_community)
            community_validity=validateCommunity(starter_community,args.edit_value,dict_ed)
            if community_validity:
                print("Community is not valid due to the following members:")
                it = iter(community_validity)
                for x in it:
                    print ("{},{}".format(x, next(it)))
                sys.exit()
            else:
                print("Community is valid")
                community=[]
                community=starter_community
                print("Starting community with:{}".format(community))
                community=loopforCommunity(community,args.edit_value,dict_ed)
                print("The number of community members is {}".format(len(community)))
        else:
            community=withoutcommunityinput(dict_ed,args.edit_value)
            print("Starting community with:{}".format(community))
            community=loopforCommunity(community,args.edit_value,dict_ed)
            print("The number of community members is {}".format(len(community)))
    else:
        sys.exit("ERROR:No edit distance given")
    if (args.strain):
        strain_info=straininfo(args.strain)
        if (args.output):
            print("Joining information")
            print("Writing output file")
            joinstrain(strain_info,community,args.output)
        else:
            print("Creating output file")
            outputfile=outfile(args.edit_value)
            print("Joining information")
            joinstrain(strain_info,community,outputfile)
    else:
        if (args.output):
            print("Writing output file")
            outputnostrain(args.output,community)
        else:
            print("Creating output file")
            outputfile=outfile(args.edit_value)
            print("Writing output file")
            outputnostrain(outputfile,community)

def correct(args):
    pass

def subsample(args):
    pass

def sequenceDictionary(input_alignment):
    sequence_dict={}
    for record in SeqIO.parse(input_alignment, "fasta"):
        sequence_dict[record.id]=str(record.seq)
    return sequence_dict

def TrimPrimers(input_alignment,primer_length):
    sequence_dict={}
    for record in SeqIO.parse(input_alignment, "fasta"):
        record.seq=record.seq[primer_length:]
        record.seq=record.seq[:-primer_length]
        sequence_dict[record.id]=str(record.seq)
    return sequence_dict

nucleiccodedicitonary={"N":set(["A","G","C","T","U"]),
                    "R":set(["A","G"]),
                    "Y":set(["C","T","U"]),
                    "S":set(["C","G"]),
                    "W":set(["A","T","U"]),
                    "K":set(["G","T","U"]),
                    "M":set(["A","C"]),
                    "B":set(["C","G","T","U"]),
                    "V":set(["A","C","G"]),
                    "D":set(["A","G","T","U"]),
                    "H":set(["A","C","T","U"])}

standardcode=["A","C","G","T","U","-"]

def customeditdistance(seq1,seq2):
    edvalue=0
    for char1, char2 in zip(seq1,seq2):
            if char1==char2:
                edvalue+=0
            else:
                if char1 in standardcode:
                    if char2 in standardcode:
                        edvalue+=1
                    else:
                        if char1 in nucleiccodedicitonary[char2]:
                            edvalue+=0
                        else:
                            edvalue+=1
                else:
                    if char2 in standardcode:
                        if char2 in nucleiccodedicitonary[char1]:
                            edvalue+=0
                        else:
                            edvalue+=1
                    else:
                        if len(nucleiccodedicitonary[char1].intersection(nucleiccodedicitonary[char2])):
                            edvalue+=0
                        else:
                            edvalue+=1
    return edvalue

def editDistanceDictionary(sequence_dict):
    dict_ed={}
    for key,value in sequence_dict.items():
        for key1,value1 in sequence_dict.items():
            if key == key1:
                continue
            else:
                edvalue=customeditdistance(value,value1)
                if (key in dict_ed.keys()):
                    if (edvalue in dict_ed[key].keys()):
                        dict_ed[key][edvalue].append(key1)
                    else:
                        dict_ed[key][edvalue]=[key1]
                else:
                    dict_ed[key]={}
                    dict_ed[key][edvalue]=[key1]
    return dict_ed


def startcommunity(input_community_file):
    starter_community=[]
    for line in input_community:
        line=line.strip()
        starter_community.append(line)
    return starter_community

def validateCommunity(input_community_file,editdistance_value,edit_distance_dictionary):
    community_validity=[]
    distances=list(range(editdistance_value)) #Create a list of numbers below edit distance value
    distances.append(editdistance_value)# append edit distance value to list
    for pair in combinations(input_community_file,2):
        if any([pair[0] in edit_distance_dictionary[pair[1]][editdistance_value] for editdistance_value in distances if editdistance_value in edit_distance_dictionary[pair[1]]]):
            community_validity.append(pair[0])
            community_validity.append(pair[1])
        else:
            continue
    return community_validity

def withoutcommunityinput(edit_distance_dictionary,editdistance_value):
    community = []
    smallest = 500 #number used that is bigger than number of possibilities
    EDnot_in_dict=[]# empty list but will contain members with no values at edit distance
    for key in edit_distance_dictionary: #loop through keys in the dictionaty made above
        if editdistance_value not in edit_distance_dictionary[key]:#Check if the key has an editdistance of input
            EDnot_in_dict.append(key) #if not append the list
        else: # if it does have that edit distqance
            if len(edit_distance_dictionary[key][editdistance_value])<smallest:#Check if the length of the values for the key at edit distance of input is less than the smallest number
                smallest= len(edit_distance_dictionary[key][editdistance_value])#set smallest to the length of the smallest list
                smallest_sequence = key#Set smallest seqeunce equal to the key
    if not EDnot_in_dict:# check if this list is empty
        community.append(smallest_sequence) # if it is empty append community with smallest sequence
    else: # if not empty
        community.append(random.choice(EDnot_in_dict))# choose random sequence from list to append community
    print(community)
    print(len(EDnot_in_dict))
    return community

def loopforCommunity(community,editdistance_value,edit_distance_dictionary):
    not_community=[]#list of members to not put in community
    distances=list(range(editdistance_value)) #Create a list of numbers below edit distance value
    distances.append(editdistance_value)# append edit distance value to list
    while True:# while there are memebers to loop through
        smallest=500
        smallest_sequence=""
        EDnot_in_dict=[]
        for key in edit_distance_dictionary:# loop through keys agains
            member = False # set to false to keep looping
            if key in community:#check if key is already in the communitiy
                continue #if key in community skip
            else:#if key is not in community check if it is the value for the community members
                for key2 in community:# loop through keys in community
                    if any([key in edit_distance_dictionary[key2][editdistance_value] for editdistance_value in distances if editdistance_value in edit_distance_dictionary[key2]]):
                        #check if key you are examining is a member of any member of the communities values at the edit distance 0 to input
                        member=True # set member to true
                        if not key in not_community: # if that key is not in not_community list add it
                            not_community.append(key)
                        break# if it is true break this loop
                if not member: #if key is not a member
                    if editdistance_value in edit_distance_dictionary[key]:# check if edit value is  present for key
                        if len(edit_distance_dictionary[key][editdistance_value])<smallest:#Check if the length of the values for the key at edit distance of input is less than the smallest number
                            smallest= len(edit_distance_dictionary[key][editdistance_value])#set smallest to the length of the smallest list
                            smallest_sequence= key#Set smallest seqeunce equal to the key
                    else: # if edit value isnt present in key
                        EDnot_in_dict.append(key) # add to list
        if smallest_sequence or EDnot_in_dict: # if either of these have a value
            if EDnot_in_dict: # if this list is not empty
                community.append(random.choice(EDnot_in_dict)) # choose random member from list
            else: # if it is empty
                community.append(smallest_sequence) # append community with smallest sequence
        else: # if they dont have a value
            break # break because we are out of members

    return community

def straininfo(metadata):
    strain_info ={} # empty dicionary for strain information
    for i, line in enumerate(metadata): # count lines in file
        line=line.strip()
        if (not i == 0): # if line is not the 1st (so the header)
            fields=line.split("\t") # split on tabs
            strain_info[fields[-1]]=fields[0]
    return strain_info

def outfile(editdistance_value):
    outputfile=open('Community_ED{}.txt'.format(editdistance_value),"w+")
    return outputfile

def joinstrain(metadata,community,outfile_name):
    for key, value in metadata.items(): # loop through keys and values in strain dictionary
        if key in community: # if that key is found in community
            outfile_name.write("{}\t{}\n".format(key, value))

def outputnostrain(outfile_name,community):
    for i in community:
        outfile_name.write("{}\n".format(i))

def outputfasta():
    pass

if __name__ == "__main__":
    main()
