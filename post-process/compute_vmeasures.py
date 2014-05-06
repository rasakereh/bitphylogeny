import os
import sys
import csv
import numpy 
from sklearn import metrics
from util import *


def get_vmeasure_traces(filepath, true_label):
    fn = [x for x in os.listdir(filepath) if 'label_traces.csv' == x][0]
    reader = csv.reader(open(filepath+'/'+fn,'rb'), delimiter=',')
    label_traces = numpy.array( list(reader) ).astype('float').astype('int')
    vmeasure_traces = numpy.array([metrics.homogeneity_completeness_v_measure
                               (label, true_label) for label in label_traces])
    header = ['homogeneity','completeness', 'v_measure']
    write_traces2csv(filepath+'/vmeasure_traces.csv',vmeasure_traces, header)

def get_true_label(fn):
    reader = csv.reader(open(fn,'rb'), delimiter=',')
    true_label = numpy.array( list(reader) )[1:,8].astype('float').astype('int')
    return (true_label)

def load_best_label(filepath):
    fn = [x for x in os.listdir(filepath) if 'best_label.csv' == x][0]
    reader = csv.reader(open(filepath+'/'+fn,'rb'), delimiter=',')
    best_label = numpy.array( list(reader) )[1:].astype('float').astype('int')
    return(best_label)


## compute mpear label v-measure

path = sys.argv[1]
traces = [x for x in os.listdir(path) if 'mcmc-traces' == x][0]
fp = path+'/'+traces
best_label = load_best_label(fp)
true_label = get_true_label( sys.argv[2] )
vmeasure = metrics.homogeneity_completeness_v_measure(best_label[:,0], true_label)
vmeasure = numpy.array(vmeasure)[numpy.newaxis,:]
header = ['homogeneity','completeness', 'v_measure']
write_traces2csv(fp + '/'+'best_label_vmeasure.csv',vmeasure, header)

## compute v-measure traces
get_vmeasure_traces(fp, true_label)
