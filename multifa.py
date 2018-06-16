from flask import Flask, request, send_file, make_response
import logging

app = Flask(__name__)

virtual_fai = None
from intervaltree import IntervalTree
virtual_fai_offsets = IntervalTree()

@app.route('/virtualgenome.fa.fai')
def return_fai():
    return get_fai()

@app.route('/virtualgenome.fa')
def return_fa_slice():
    import sys
    for header in request.headers.keys():
        sys.stderr.write(header + " " + request.headers.get(header) + "\n")

    byterange = request.headers.get('Range')
    sys.stderr.write("byterange="+str(byterange))
    return get_fa_slice(byterange)

def get_fai():
    """constuct a virtual genome by aggregating all the individual genomic fais
    if there is none yet for the current application context.
    """
    global virtual_fai
    if not virtual_fai:
        virtual_fai = create_virtual_fai()
    return virtual_fai

def create_virtual_fai():
    """implementation details of virtual fai construction
    """
    retval = ''
    offset = 0
    global virtual_fai_offsets
    #for now, just to get something working
    for gnm in ('data/mygenome1.fa.fai', 'data/mygenome2.fa.fai'):
        thisfile = open(gnm)
        for line in thisfile: 
            data = line.split('\t')
            new_offset = int(data[2]) + offset
            data[2] = str(new_offset)
            retval += '\t'.join(data)
        thisfile.close()
        #figure out offset (ie what the concatenated offset would be in the virtual genome)
        lastseq_offset = int(data[2])
        lastseq_len = int(data[1])
        lastseq_line_len = int(data[3])
        lastseq_linedelim_len = int(data[4])-int(data[3])
        import math
        new_offset = offset + lastseq_offset + lastseq_len + lastseq_linedelim_len*(int(math.ceil(float(lastseq_len)/lastseq_line_len)))
        virtual_fai_offsets[offset:new_offset] = {'file':gnm[:-4], 'offset':offset}
        offset = new_offset
    return retval

def get_fa_slice(byterange):
    import re
    prog = re.compile('bytes=(\d+)-(\d+)')
    m = prog.match(byterange)
    start = int(m.group(1))
    stop = int(m.group(2))
    d = virtual_fai_offsets[start:stop].pop().data
    f = open(d['file'])
    f.seek(start-d['offset'])
    r = make_response(send_file(f, 'X-application-fasta'), 206)
    r.headers['Content-Length'] = stop-start+1
    r.headers['Accept-Ranges'] = 'bytes'
    return r

