import os
import csv
import numpy 
from sklearn import metrics

def write_traces2csv(file_name, traces, header=None):
    fh = open(file_name, 'wb')
    writer = csv.writer(fh)
    if header:
        writer.writerow(header)
    [writer.writerow(x) for x in traces]
    fh.close()

def get_vmeasure_traces(filepath, true_label):
    fn = [x for x in os.listdir(filepath) if 'label_traces' in x][0]
    reader = csv.reader(open(filepath+'/'+fn,'rb'), delimiter=',')
    label_traces = numpy.array( list(reader) ).astype('float').astype('int')
    vmeasure_traces = numpy.array([metrics.homogeneity_completeness_v_measure
                               (label, true_label) for label in label_traces])
    header = ['homogeneity','completeness', 'v_measure']
    write_traces2csv(filepath+'/v_measure_traces.csv',vmeasure_traces, header)


def find_synth_folders(pathname):
    codename_folders = os.listdir(pathname)
    filepath=[]
    for f in codename_folders:
        samplefolder = os.listdir(pathname+'/'+f)[0] 
        if 'noisy' in samplefolder:
            filepath.append( pathname +'/'+ f+'/'+samplefolder )
    return(filepath)


def get_true_label(filepath):
    fn = [x for x in os.listdir(filepath) if '0.01' in x][0]
    reader = csv.reader(open(filepath+'/'+fn,'rb'), delimiter=',')
    true_label = numpy.array( list(reader) )[1:,8].astype('float').astype('int')
    return (true_label)

pathname = '/home/yuan03/Downloads/traces-eth/Sottoriva2/pairwise/'
filepath = find_synth_folders(pathname)
datafolder = '/home/yuan03/git-repo/phylo-tree/code/data/full_methy/'
for fp in filepath:
    print 'Processing ', fp
    true_label = get_true_label(datafolder)
    get_vmeasure_traces(fp, true_label)

 
