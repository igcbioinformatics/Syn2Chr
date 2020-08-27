#!/usr/bin/python3
# Version 1(2018-7-30):
# Build pseudo-chromosome from scaffolds, downstream of SyntenyBuild.
# chrbuild ScaffoldInput.fa chr [-]scaffoldname, [-]scaffoldname, ...
# '-' means complementary, '+' means gaps. K/M/G accepted.
# Example: chrbuild.py mOncTor1.PB.asm4.scaff2.break1.gap1.polish2.fa Chr11 +3M,Super-Scaffold_3029,+100K,Super-Scaffold_6058,+100K,Super-Scaffold_334,+10K,-Super-Scaffold_444,Contig161arrow_subseq_1:222120_obj,+50K,Super-Scaffold_622,+50K,Super-Scaffold_3028 > Chr11_mOncTor_asm4.2112.fa
# Version 2(2018-10-2):
# Add function of Gold Path .AGP 2.0 output 
# Version 2(2018-10-4):
# Add index checking, add AGP head suppressor
# Version 2019-3:
# Add AGP only output, add 'AGP'to the begining of Synteny will supress sequence output
# Version 2019-5:
# Add arguments by introducing the argparse module. Prepare script for replacing the long command line input by using input synteny file.
# Version 2019-6:
# Introduce the ability to create AGP files, fasta files and summary files.
# Version 2019-7:
# Restructure code to improve performance, add option to split scaffolds/contigs. 
# Jl17 optimize
# Version 2020-4:
# Add '@' option to synteny file
import os
import sys
import argparse
from datetime import date

################################ Argparse module #############################

parser = argparse.ArgumentParser(prog='ChrBuild', description='Build pseudoc-chromosomes from synteny file. Requires python 3.4+')
parser.add_argument('Scaffold_FASTA', help = 'Input file name Scaffold fasta file.')
parser.add_argument('Synteny_File', help = 'Input synteny file for Pseudo chromosome building. Synetny file ')
parser_chr = parser.add_argument_group('ChrBuilder - functions:')
parser_chr.add_argument('-S','--Summary', dest='summary', action = 'store_true', default = False, help = 'Instead of all three (= default), only the summary file is created')
parser_chr.add_argument('-A','--AGP', dest='agp', action = 'store_true', default = False, help = 'Instead of all three (= default), only the AGP file is created')
parser_chr.add_argument('-F','--Fasta', dest='sequence', action = 'store_true', default = False, help = 'Instead of all three (= default), only the fasta file is created')
parser_chr.add_argument('-L','--Liftover', dest='gff', action = 'store', nargs = '?', const = True, default = '', help = 'Liftover gtf/gff3 annotation from scaffolds to genome based on synteny file.')
parser_chr = parser.add_argument_group('ChrBuilder - settings:')
parser_chr.add_argument('-D', dest='split_symbol', action = 'store', default = ':', metavar='SYMBOL', help = 'Change the contig split symbol. Default split symbol ":".')
parser_chr.add_argument('-H','--Nohead', action = "store_true", default = False, help = "Remove header from summary and AGP file.")
parser_chr.add_argument('-o', dest= 'Output', action = 'store', metavar='NAME', nargs='?', const = '', default = '', help = 'Give a custom name to the output file(s). Do not use suffixes.')
parser_chr.add_argument('-t', dest= 'Toplevel', action = 'store_true', default = False, help = 'Output Toplevel sequence of fasta. All unused scaffolds will be outputed.')
args = parser.parse_args()

################################## script ####################################
# Initiate:
agp_requested, sequence_requested, summary_requested, header_set, split_symbol, top_level = args.agp, args.sequence, args.summary, not args.Nohead, args.split_symbol, args.Toplevel
Chr_list, summary_total, Used_scaf, summary_output = [],0,[],'\t'.join(['Chromosome','Scaffold_length','Scaffold_count'])+'\n'
filename = args.Output if args.Output != '' else os.getcwd().split('/')[-1]
if not (agp_requested or sequence_requested or summary_requested): agp_requested = sequence_requested = summary_requested = True
if args.gff == True:
    print('Option Liftover is given, but does not following by gtf/gff/gff3 file name. Liftover ommited.' )
    args.gff = '' 
# Find and read files
Fas = open(args.Scaffold_FASTA,'r')
try: Faindex = open(args.Scaffold_FASTA+'.fai','r')
except FileNotFoundError:
    sys.stderr.write('ERROR: FASTA index file '+ args.Scaffold_FASTA +'.fai can not be found.')
    exit(1)
Fai = {dx.split('\t')[0]:list(map(int,dx.split('\t')[1:])) for dx in Faindex} 
Faindex.close()
with open(args.Synteny_File,'r') as SCL: Scaf_list = [dx.strip() for dx in SCL.readlines() if dx != '\n' and dx.strip()[0] != '#']

#Header info
header_info = Scaf_list[0].split('\t')
if header_info[3] == '-': header_info[3] = date.today().strftime('%d-%B-%y')
agp_output = '##agp-version\t2.0\n#ORGANISM:\t'+header_info[0]+'\n#TAX_ID:\t'+header_info[1]+'\n#ASSEMBLY NAME:\t'+header_info[2]+'\n#ASSEMBLY DATE:\t'+header_info[3]+'\n#GENOME CENTER:\t'+header_info[4]+'\n#DESCRIPTION:\t'+header_info[5]+'\n#COMMENTS:\t'+header_info[6]+'\n' if header_set else ''
for dx in Scaf_list[1:]:
    if ' ' not in dx : print(dx, '\ndoes not follow synteny line syntex. Deleted'); continue
    Chr_list.append(dx.split(' '))
# Main:
if sequence_requested: output = open(filename+'.fa','w')
for dx in Chr_list:
    Chr_length, summary_length,scaffold_count, Chr_output, line_count, Chr_pointer = 0,0,0,'', 1, 1
    print('Processing {}'.format(dx[0]))
    for scaffoldname in dx[1].split(','):
        # Add gaps when requested
        if scaffoldname[0] in '+@':
            seq_length = round(int(scaffoldname[1:-1])*{'G':1E9,'M':1E6,'K':1E3}[scaffoldname[-1].upper()]) if scaffoldname[-1] in 'GMKgmk' else round(int(scaffoldname[1:]))
            if scaffoldname[0] == '@': 
                if Chr_length > seq_length: print('WARNING: @{} will overlap with previous contig(s).'.format(scaffoldname[1:]))
                seq_length = seq_length - Chr_length - 1
            agp_output += '\t'.join([dx[0], str(int(Chr_pointer)),str(int(Chr_pointer+seq_length - 1)), str(line_count),'N',str(seq_length),'scaffold','no','na\n'])
            if sequence_requested: Chr_output += 'N'*seq_length
            Chr_length += seq_length
        else:
            ## Check for reverse complementary scaffolds
            if scaffoldname[0] == '-':
                rc_flag = True
                scaffoldname = scaffoldname[1:]
            else: rc_flag = False
            ## Check for presence of broken scaffold symbol:
            if split_symbol in scaffoldname:
                name_split = scaffoldname.split(split_symbol)
                if len(name_split) != 2 or len(name_split[1].split('-')) != 2:
                    sys.stderr.write('ERROR: The scaffold {} containing scaffold split symbol "{}". But the format is incomprehensible(Scaffld:start-end expected, e.g. Chr1:1-100K). Use -D option to change default split symbol if necessary'.format(scaffoldname, split_symbol))
                    exit(1)
                scaffoldname = name_split[0]
                split_start, split_end = [int(cx[:-1])*{'K':1e3,'M':1e6,'G':1e9}[cx[-1].upper()] if cx[-1] in 'GMKgmk' else int(cx) for cx in name_split[1].split('-')]
                if split_start > split_end: 
                    sys.stderr.write('WARNING: Given region in scaffold {} looks reversed. Please read the instruction on spliting reverse complementary scaffolds.'.format(scaffoldname))
                    split_start, split_end = split_end, split_start
                broken_contig_flag = True
            else:
                broken_contig_flag = False

            ## check if fai file is correct
            if not broken_contig_flag and scaffoldname not in Fai.keys():
                sys.stderr.write('ERROR: Following scaffold is not found in {}: {}, please check if fai index file is correct.\n'.format(dx[0],scaffoldname))
                exit(1)
            if scaffoldname in Used_scaf: sys.stderr.write('WARNING: {} used more than one times in the synteny.\n'.format(scaffoldname))
            if sequence_requested: Fas.seek(Fai[scaffoldname][1])
            if rc_flag:
                seq_scaffold = Fas.read(Fai[scaffoldname][0] + (Fai[scaffoldname][3]-Fai[scaffoldname][2]) * (Fai[scaffoldname][0] // Fai[scaffoldname][2])).replace('\n','').lower()[::-1]
                for ex,cx in [['a','T'],['t','A'],['g','C'],['c','G'],['n','N']]:
                    seq_scaffold = seq_scaffold.replace(ex,cx)
            else:
                seq_scaffold = Fas.read(Fai[scaffoldname][0] + (Fai[scaffoldname][3]-Fai[scaffoldname][2]) * (Fai[scaffoldname][0] // Fai[scaffoldname][2])).replace('\n','').upper()
            if broken_contig_flag:
                seq_scaffold = seq_scaffold[split_start:split_end]
            Chr_output += seq_scaffold
            Used_scaf.append(scaffoldname)
            seq_length = split_end-split_start + 1 if broken_contig_flag else Fai[scaffoldname][0]
            orientation = '-\n' if rc_flag else '+\n'
            name = scaffoldname + '[' + str(split_start) + '-' + str(split_end) + ']' if broken_contig_flag else scaffoldname
            agp_output += '\t'.join([dx[0], str(int(Chr_pointer)),str(int(Chr_pointer+seq_length - 1)), str(line_count),'W', name,'1',str(seq_length),orientation])
            scaffold_count += 1
            summary_length += seq_length
            Chr_length += seq_length
        Chr_pointer += seq_length
        line_count += 1
    if sequence_requested:#slice of string can  be out of range, return empty. this solves last line shorter than 60bp
        output.write('>' + dx[0] + '\n' + '\n'.join(Chr_output[i:i+60] for i in range(0, len(Chr_output), 60)) + '\n')
    if summary_requested:# summary section
        summary_output += '\t'.join([dx[0], str(Chr_length),str(summary_length),str(scaffold_count)])+'\n'
        summary_total += summary_length
if top_level:
    print('Processing unplaced scaffolds.')
    for dx in Fai.keys():
        if dx not in Used_scaf:
            agp_output += '\t'.join([dx,'1',str(Fai[dx][0]),'1','W', dx,'1',str(Fai[dx][0]),'+'])+'\n'
            if sequence_requested:
                Fas.seek(Fai[dx][1])
                seq_scaffold = Fas.read(Fai[dx][0] + (Fai[dx][3]-Fai[dx][2]) * (Fai[dx][0] // Fai[dx][2])).replace('\n','')
                output.write('>' + dx + '\n' + '\n'.join(seq_scaffold[cx:cx+60] for cx in range(0, len(seq_scaffold), 60)) + '\n')
Fas.close()
if sequence_requested: output.close()
if agp_requested:
    with open(filename+'.agp','w') as output: output.write(agp_output)
if summary_requested:
    with open(filename+'.txt','w') as output:
        output.write(summary_output)
        scaffold_total = sum(Fai[dx][0] for dx in Fai.keys())
        output.write('\t'.join(['Assembled:',str(summary_total),str(len(Used_scaf))])+'\n')
        output.write('\t'.join(['Total:',str(scaffold_total),str(len(Fai))])+'\n')
        output.write('Sequence Assembled: {0:3.1f}%'.format(summary_total/scaffold_total*100))

# Function gff3 liftover
if args.gff != '':
    AGP_from, AGP_to, gff3, used_new = {}, {}, [], []
    #AGP_from[scaffold_name] = [start, end, scaffold_name, direction]
    AGP_from = dict(zip(Fai.keys(),[[1, Fai[dx][0],dx,'+'] for dx in Fai.keys()]))
    for line in agp_output:
        if line.strip()[0] == '#' or '\t' not in line: continue
        dx = line.strip().split('\t')
        if dx[0] not in used_new : used_new.append(dx[5])
        AGP_to[dx[5]] = [dx[0],int(dx[1]),int(dx[2]),dx[8]]
    for line in open(args.gff,'r'):
        if line.strip()[0] == '#' or '\t' not in line:
            print(line)
            continue
        dx = line.rstrip().split('\t')
        dx[3],dx[4] = int(dx[3]), int(dx[4])
        if dx[0] not in used_new and dx[0] not in AGP_from.keys():
            gff3.append(dx)
            continue
        if dx[0] in used_new and dx[0] not in AGP_from.keys():
            if AGP_to[dx[0]][3] == '+':
                pos1 = AGP_to[dx[0]][1] + dx[3] -1
                pos2 = AGP_to[dx[0]][1] + dx[4] -1
                dir_ = dx[6]
            else:
                pos1 = AGP_to[dx[0]][2] - dx[4] +1
                pos2 = AGP_to[dx[0]][2] - dx[3] +1
                dir_ = {'+':'-','-':'+'}[dx[6]]
            gff3.append([AGP_to[dx[0]][0],dx[1],dx[2],pos1,pos2,dx[5],dir_,dx[7],dx[8]])
            print(AGP_to[dx[0]][0],dx[1],dx[2],pos1,pos2,dx[5],dir_,dx[7],dx[8], file = sys.stderr)
            continue
        hit1, hit2 = False, False
        chr_1, chr_2 = 'na','na'
        for start_,end_,scaffold,direct in AGP_from[dx[0]]:
            if direct == 'na': continue
            if start_ <= dx[3] and dx[3] <= end_:
                if direct == AGP_to[scaffold][3]:
                    pos1 = dx[3] - start_ + AGP_to[scaffold][1]
                    dir_ = dx[6]
                else:
                    pos2 = AGP_to[scaffold][2] - dx[3] + start_
                    dir_ = '|' if dx[6] == '+' else 'x'
                chr_1 = AGP_to[scaffold][0]
                hit1 = True
            if start_ <= dx[4] and dx[4] <= end_:
                if direct == AGP_to[scaffold][3]:
                    pos2 = dx[4] - start_ + AGP_to[scaffold][1]
                    dir_ = dx[6]
                else:
                    pos1 = AGP_to[scaffold][2] - dx[4] + start_
                    dir_ = '|' if dx[6] == '+' else 'x'
                chr_2 = AGP_to[scaffold][0]
                hit2 = True
            if hit1 and hit2 : break
        if hit1 and chr_2 =='na':
            hit2 = hit1 -dx[3] + dx[4] if dir_ in ['+','-'] else hit1 +dx[3] - dx[4] -1
            chr_2 = chr_1
            print('Warning:',dx[8],'on',dx[0],"have 3' end fall inside scaffold gaps, please check", file=sys.stderr)
        if hit2 and chr_1 =='na':
            hit1 = hit2 -dx[4] + dx[3] if dir_ in ['+','-'] else hit2 -dx[3] + dx[4] -1
            chr_1 = chr_2
            print('Warning:',dx[8],'on',dx[0],"have 5' end fall inside scaffold gaps, please check", file=sys.stderr)
        if hit1 and hit2 and chr_1 == chr_2:
            gff3.append([chr_1,dx[1],dx[2],pos1,pos2,dx[5],dir_,dx[7],dx[8]])
        else:
            print(dx, 'error')
            exit(1)
    # direction fix and output:
    exon_count, CDS_count, gene = 0,0,[]
    pointer = 0
    line_count = len(gff3)
    while True:
        dx = gff3[pointer]
        if dx[2] == 'gene':
            print('\t'.join(map(str,[dx[0],dx[1],dx[2],dx[3],dx[4],dx[5],{'+':'+','-':'-','x':'+','|':'-'}[dx[6]],dx[7],dx[8]])))
            pointer += 1
            continue
        elif dx[2] == 'mRNA':
            print('\t'.join(map(str,[dx[0],dx[1],dx[2],dx[3],dx[4],dx[5],{'+':'+','-':'-','x':'+','|':'-'}[dx[6]],dx[7],dx[8]])))
            pointer += 1
            gene_name = dx[8].split('ID=')[1].split(';')[0]
            continue
        if dx[2] == 'exon': exon_count += 1
        if dx[2] == 'CDS': CDS_count += 1
        gene.append(dx)
        pointer +=1
        if gff3[pointer][2] in ['gene','mRNA'] or pointer+1 >= line_count:
            if dx[6] in ['|','x']:
                for cx in range(len(gene)):
                    if gene[cx][2] == 'exon':
                        gene[cx][8] = 'ID'+gene_name+'.exon'+ str(exon_count-int(gene[cx][8].split('.exon')[1].split(';')[0])+1) + ';Parent=' + gene_name
                    elif gene[cx][2] == 'CDS':
                        gene[cx][8] = 'ID'+gene_name+'.CDS'+ str(exon_count-int(gene[cx][8].split('.CDS')[1].split(';')[0])+1) + ';Parent=' + gene_name
            gene = sorted(gene, key = lambda ex:ex[3])
            gene = sorted(gene, key = lambda ex:ex[4])
            for cx in gene:
                print('\t'.join(map(str,[cx[0],cx[1],cx[2],cx[3],cx[4],cx[5],{'+':'+','-':'-','x':'+','|':'-'}[cx[6]],cx[7],cx[8]])))
            exon_count, CDS_count, gene = 0,0,[]
            #print('Processing:',gene_name, file=sys.stderr)
            if pointer+1 >= line_count: exit(0)
print('\nDONE')
exit(0)
