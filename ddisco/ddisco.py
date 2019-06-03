import sys
import re
import random
import os
import argparse
import functools
import operator
import copy
from itertools import combinations
from collections import defaultdict
from Bio import SeqIO
from ddisco._version import __version__

def main():
    parser = argparse.ArgumentParser(prog="ddisco", add_help=False)
    parser.add_argument('-v', '--version', action='version', version="v{}".format(__version__))
    parser.add_argument("-i", "--input", type=argparse.FileType("r"),dest="input",
                      help="alignment file in fasta form")
    parser.add_argument("--editdistance",type=int, dest="edit_value",
                     help="Edit distance value as interger")
    parser.add_argument("-o","--output",type=argparse.FileType("w"),dest="output",
                        help="output file name")
    parser.add_argument("--seed",type=int,default=os.urandom(64),
                        help="seed number for reproducibility")


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

    #Subsample subcommand
    parser_subsample=subparsers.add_parser("subsample", parents=[parser],help="Module to subsample highly diverse community")
    parser_subsample.add_argument("--num-taxa", "-n", type=int, dest="num_taxa",
                    help="number of strains desired in final community")
    parser_subsample.add_argument("--community", "-c", type=argparse.FileType("r"), 
                                  help="Tab seperated file with taxa ids in the first column with metadata in additional columns")
    parser_subsample.add_argument("--group-by", "-g", dest="group_by",
                        help="Column name to group-by for proportion calculation. Default to second column")
    parser_subsample.add_argument("--proportion", "-p", type=argparse.FileType("r"),dest="proportion",
                        help="File of the relative proportions of each taxonomic rank desired in final community")
    parser_subsample.add_argument("--taxa-num-enforce", "-e", action="store_true", dest="num_enforce", default=False, help="")
    parser_subsample.set_defaults(func=subsample)

    # Parse args
    args = parser.parse_args()
    args.func(args)

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

def subsample(args):
    #parser_subsample=subparsers.add_parser("subsample", parents=[parser],help="Module to subsample highly diverse community")
    #parser_subsample.add_argument("--num-taxa", "-n", type=int, dest="num_taxa",
    #                help="number of strains desired in final community")
    #parser_subsample.add_argument("--community", "-c", type=argparse.FileType("r"), 
    #                              help="Tab seperated file with taxa ids in the first column with metadata in additional columns")
    #parser_subsample.add_argument("--group-by", "-g", dest="group_by",
    #                    help="Column name to group-by for proportion calculation. Default to second column")
    #parser_subsample.add_argument("--proportion", "-p", type=argparse.FileType("r"),dest="proportion",
    #                    help="File of the relative proportions of each taxonomic rank desired in final community")
    #parser_subsample.add_argument("--taxa-num-enforce", "-e", action="store_true", default=False, help="")
    if args.num_enforce and not args.proportion:
        print("ERROR: The --taxa-num-enforce option must be used with the --proportion option", file=sys.stderr)
        sys.exit()
    elif args.num_enforce and not args.num_taxa:
        print("ERROR: The --taxa-num-enforce option must be used with the --num-taxa option", file=sys.stderr)
        sys.exit()
    elif not (args.num_taxa or args.proportion):
        print("ERROR: Either the --num-taxa or --proportion option required", file=sys.stderr)
        sys.exit()


    random.seed(args.seed)
    community_list = [line.strip().split("\t") for line in args.community]
    header = community_list[0]
    community_list = community_list[1:]
    #Check to see num_taxa great or equal to current community size
    if args.num_taxa and len(community_list) <= args.num_taxa:
        print("ERROR: --num-taxa >= current community taxa count", file=sys.stderr)
        sys.exit()

    if args.num_taxa and not args.proportion:
        new_community=random.sample(community_list,args.num_taxa)
        with open("Subsampled_community_taxa{}.txt".format(args.num_taxa),"w") as subsample_output:
            print("\t".join(header), file=subsample_output)
            for member in new_community:
                print("{}".format("\t".join(member)), file=subsample_output)
    elif args.proportion:
        group_by_col = 1
        if args.group_by:
            if args.group_by in header:
                group_by_col = header.index(args.group_by)
            else:
                print("ERROR: Group by variable ({}) was not found in the header. The provided headings are {}".format(group_by, ",".join(header)), file=sys.stderr)
                sys.exit()

        grouping_dict = defaultdict(list)
        for taxon in community_list:
            grouping_dict[taxon[group_by_col]].append(taxon)
        
        goal_prop = defaultdict(float)
        for prop in args.proportion: 
            prop = prop.strip()
            prop_split = prop.split("\t")
            goal_prop[prop_split[0]] = float(prop_split[1])
            if not prop_split[0] in grouping_dict:
                print("ERROR: {} not found in grouping column".format(prop[0]), file=sys.stderr)
                sys.exit()
                
        if not isclose(1, sum(goal_prop.values())):
            print("ERROR: Proportions do not sum to 1. Currently sum to {}".format(sum(goal_prop.values())), file=sys.stderr)
            sys.exit()

        grouping_dict = {k: grouping_dict[k] for k in goal_prop}
        #shuffle taxa
        for group in grouping_dict:
            random.shuffle(grouping_dict[group])
        current_total = sum(len(lst) for lst in grouping_dict.values())
        current_props = {k: len(grouping_dict[k])/current_total for k in grouping_dict}
        current_sse = 0
        error_dict = defaultdict(float)
        for group in current_props:
            error_dict[group] = current_props[group] - goal_prop[group]
            current_sse += (current_props[group] - goal_prop[group])**2

        #loop to minimize sum of squared "errors"
        while True:
            #find group with maximum difference to goal proportion
            max_error = max(error_dict.items(), key=operator.itemgetter(1))
            #if max diffrence is zero or negative can not improve and break
            #or if max_error group only has one member.
            if max_error[1] <= 0 or len(grouping_dict[max_error[0]]) == 1:
                break
            #if args.num_enforce check to make sure removing doesn't drop taxa below args.num_taxa
            if args.num_enforce and (current_total - 1) < args.num_taxa:
                break
            #Remove taxa from group with max error
            test_grouping_dict = {}
            for group in grouping_dict:
                if group == max_error[0]:
                    temp_list = copy.deepcopy(grouping_dict[group])
                    temp_list.pop()
                    random.shuffle(temp_list)
                    test_grouping_dict[group] = copy.deepcopy(temp_list)
                else:
                    test_grouping_dict[group] = copy.deepcopy(grouping_dict[group])
            #Recalculate SSE
            temp_total = current_total - 1
            temp_props = {k: len(test_grouping_dict[k])/temp_total for k in test_grouping_dict}
            temp_sse = 0
            temp_error_dict = defaultdict(list)
            for group in temp_props:
                temp_error_dict[group] = temp_props[group] - goal_prop[group]
                temp_sse += (temp_props[group] - goal_prop[group])**2
            #Test if new SSE is less then current and proceed accordingly
            if temp_sse < current_sse:
                current_total -= 1
                current_props = temp_props.copy()
                current_sse = temp_sse
                error_dict = temp_error_dict.copy()
                grouping_dict[max_error[0]] = copy.deepcopy(test_grouping_dict[max_error[0]])
            else:
                break
            
        if args.num_taxa and sum(len(lst) for lst in grouping_dict.values()) > args.num_taxa:
            down_prop = 1 - (args.num_taxa/sum(len(lst) for lst in grouping_dict.values()))
            for group in grouping_dict:
                group_remove_num = round(len(grouping_dict[group])*down_prop)
                #check if num remove will drop below num_taxa if num_enforce and adjust if needed
                if args.num_enforce and (sum(len(lst) for lst in grouping_dict.values()) - group_remove_num) < args.num_taxa:
                    group_remove_num -= args.num_taxa - (sum(len(lst) for lst in grouping_dict.values()) - group_remove_num)
                if group_remove_num < 1 or len(grouping_dict[group]) == 1:
                    continue
                for _ in range(group_remove_num):
                    if len(grouping_dict[group]) == 1:
                        break
                    else:
                        grouping_dict[group].pop()
                        current_total -= 1
                if sum(len(lst) for lst in grouping_dict.values()) == args.num_taxa and args.num_enforce:
                    break
            #final correction to make sure taxa <= num_taxa if not args.num_enforce
            if sum(len(lst) for lst in grouping_dict.values()) > args.num_taxa:
                diff = sum(len(lst) for lst in grouping_dict.values()) - args.num_taxa
                while (diff > 0):
                    group_counts = {k: len(grouping_dict[k]) for k in grouping_dict}
                    group = sorted(grouping_counts.items(), key=operator.itemgetter(1), reverse=True)[0]
                    if group[1] == 1:
                        break
                    else:
                        grouping_dict[group[0]].pop()
                        current_total -= 1
                        diff -= 1



        #Recalculate proportions after final correction
        current_props = {k: len(grouping_dict[k])/current_total for k in grouping_dict}
        #print number of taxa
        print("Number of Taxa: {}".format(current_total))
        #print proportions
        print("Actualized proportions")
        for group in current_props:
            print("{}:{:0.4f}".format(group, current_props[group]), file=sys.stderr)
        with open("Subsampled_community_taxa_prop.txt", "w") as subsample_output:
            print("\t".join(header), file=subsample_output)
            for group in grouping_dict:
                for taxa in grouping_dict[group]:
                    print("{}".format("\t".join(taxa)), file=subsample_output)

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

@functools.lru_cache(maxsize=100000)
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
    dict_ed = defaultdict(lambda: defaultdict(list))
    for taxa in combinations(sequence_dict.keys(), 2):
        edvalue=customeditdistance(sequence_dict[taxa[0]], sequence_dict[taxa[1]])
        dict_ed[taxa[0]][edvalue].append(taxa[1])
        dict_ed[taxa[1]][edvalue].append(taxa[0])

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
    for sequence_id in community:
        outfile_name.write("{}\n".format(sequence_id))

def outputfasta(sequence_dict,community,editdistance_value):
    with open('Community_ED{}.fasta'.format(editdistance_value),"w+") as fasta_file:
        for sequence_id in community:
            fasta_file.write(">{}\n{}\n".format(sequence_id,sequence_dict[sequence_id]))

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

if __name__ == "__main__":
    main()