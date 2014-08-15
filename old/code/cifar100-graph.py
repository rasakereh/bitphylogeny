import os
import sys
import uuid
import cPickle

from util     import *
from tssb     import *
from logistic import *

file = sys.argv[1]
print file
fh = open(file)
tssb = cPickle.load(fh)
fh.close()
print "TSSB file loaded..."

num_images  = 10
tiling      = "2x5"
threshold   = 50
codes       = cifar100_codes(tssb.num_data)

if len(sys.argv) > 2:
    img_dir = sys.argv[2] + '/Data/CIFAR-100/images'
else:
    img_dir = os.environ['HOME'] + '/Data/CIFAR-100/images'

tmp_dir = '/tmp'
    
def montage(images):
    outname = "%s/%s.png" % (tmp_dir, uuid.uuid4())
    cmd1 = "montage -geometry 32x32+0+0 -tile %s %s %s" % (tiling, " ".join(images), outname)
    os.popen(cmd1)
    #cmd2 = "optipng -o7 %s" % (outname)
    #os.popen(cmd2)
    return outname

fh = open("%s.gdl" % (file), 'w')

def print_node(node, parent):
    num_data = node['node'].num_data()
    if num_data > threshold:
        indices   = list(node['node'].data)
        if len(indices) >= num_images:
            image_files  = map(lambda n: "%s/%05d.png" % (img_dir, indices[n]), argsort(map( lambda n: -node['node'].logprob(codes[n]), indices))[:num_images])
            montage_file = montage(image_files)
            node_name    = uuid.uuid4()

            print >>fh, """ \
node: { title:"%s" iconfile:"%s" borderstyle:invisible label:"" } \n\
""" % (node_name, montage_file)
            if parent is not None:
                print >>fh, """ \
edge: { source:"%s" target:"%s" arrowsize:8 } \n\
""" % (parent, node_name)
        else:
            node_name = parent
        for child in node['children']:
            print_node(child, node_name)
 
#print >>fh, """ \
#graph: { title:            "CIFAR 100 Graph"  \n\
#         layout_algorithm: tree               \n\
#         node.fontname:    "helvR8"           \n\
#         smanhattanedges:   yes               \n\
#         splines:           yes               \n\
#         outportsharing:    yes               \n\
#         equal_y_dist:      yes               \n\
#         yspace:            20                \n\
#         xspace:            5                 \
#"""
print >>fh, """ \
graph: { title:            "CIFAR 100 Graph"  \n\
         layout_algorithm: tree               \n\
         node.fontname:    "helvR8"           \n\
         smanhattanedges:   yes               \n\
         splines:           yes               \n\
         outportsharing:    yes               \n\
         equal_y_dist:      no                \n\
         yspace:            20                \n\
         xspace:            5                 \
"""
print_node(tssb.root, None)
print >>fh, """}"""
fh.close()
