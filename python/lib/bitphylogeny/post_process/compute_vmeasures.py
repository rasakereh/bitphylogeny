import os
import csv
import numpy 

from bitphylogeny.util import *

def write_vmeasure_traces(filepath, true_label):
    from sklearn import metrics
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

def load_mpear_label(filepath):
    fn = [x for x in os.listdir(filepath) if 'mpear_label.csv' == x]
    
    if len(fn) == 0:
        fn1 = [x for x in os.listdir(filepath) if 'label_traces.csv' == x][0]
        reader = csv.reader(open(filepath+'/'+fn1,'rb'), delimiter=',')
        label_traces = numpy.array( list(reader) ).astype('float').astype('int')
        return( cluster_with_maxpear(label_traces) ) 

    reader = csv.reader(open(filepath+'/'+fn[0],'rb'), delimiter=',')
    mpear_label = numpy.array( list(reader) )[1:].astype('float').astype('int')
    return(mpear_label[:,0])
