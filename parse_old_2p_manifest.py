# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 14:43:38 2019

@author: bnste
"""

import argparse
import os.path
import sys

delim=','

parser = argparse.ArgumentParser()
parser.add_argument('expt_info_files',nargs='+',type=argparse.FileType('r'), help='input experiment info files, in a comma-separated list')
parser.add_argument('--out_files',nargs='+', type=argparse.FileType('w'), default = [sys.stdout], help = 'output info file(s), comma-separated list')
args = parser.parse_args()

expt_info_files = args.expt_info_files
out_files = args.out_files
if (len(expt_info_files) != len(out_files)) and (len(out_files) != 1):
    raise Exception('out_files must be a list of length 1 or the same length as expt_info_files')


colnames = ['mouse','date','tif','wavesurfer','voyeur','hologram','maskdir','power','npulse','tpulse','ipi','inhaledelay','notes']
cur_line = dict(zip(colnames,['']*6))
#out_lines = [delim.join(colnames)]


nondelim = ';' if delim == ',' else ','



for reader_num,reader in enumerate(expt_info_files):
    if len(out_files) == 1:
        writer = out_files[0]
        if reader_num == 0:
            writer.write(delim.join(colnames) + '\n')
    else:
        writer = out_files[reader_num]
        writer.write(delim.join(colnames) + '\n')
    cur_line['mouse'], cur_line['date'] = os.path.splitext(os.path.basename(reader.name))[0].split('_')
    for line in reader:
        line = line.replace('\n','').strip()
     #   print(line)
        if len(line) <= 1:
           # print('starting new block')
            if cur_line[colnames[3]] != '':
                #out_lines = out_lines + [delim.join([cur_line[x] for x in colnames])] 
                writer.write(delim.join([cur_line[x] for x in colnames]) + '\n')
                writer.flush()
            for name in colnames[2:]:
                cur_line[name] = ''
            continue
        splitind = line.find(':')
        if splitind==-1:
            continue
      #  print('found good line!')
        prefix = line[:splitind].lower().strip()
        if prefix not in colnames:
            continue
      #  print('found prefix {}'.format(prefix))
        if prefix == 'notes':
            lineval = "\"" + line[splitind:].strip().replace(':',' ').replace("\"","") + "\""
        else:
            lineval = line[splitind+1:].strip().split(' ')[0].replace(delim, nondelim) #.split('\n')[0]
        cur_line[prefix.lower()] = lineval
    if cur_line[colnames[3]] != '':
        #out_lines = out_lines + [delim.join([cur_line[x] for x in colnames])] 
        writer.write(delim.join([cur_line[x] for x in colnames]) + '\n')
        writer.flush()
    for name in colnames[2:]:
        cur_line[name] = ''
    
        
[reader.close() for reader in expt_info_files];
[writer.close() for writer in out_files];

