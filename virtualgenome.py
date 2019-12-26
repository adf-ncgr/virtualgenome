from flask import Flask, request, send_file, make_response
import logging
import urllib2
from struct import pack, unpack

app = Flask(__name__)

virtual_fai = None
virtual_gzi = None
genome_sizes = {}
#for VirtualFile handoff strategy
genome_succession = {}
from intervaltree import IntervalTree
virtual_fai_offsets = IntervalTree()
virtual_gzi_offsets = IntervalTree()

import re
byterange_matcher = re.compile('bytes=(\d+)-(\d+)')
protocol_matcher = re.compile('^https?://')

from io import BufferedReader

class VirtualFile(BufferedReader):

    current_gnm = None
    current_f = None

    def __init__(self,current_gnm,first_seek,size):
        self.current_gnm=current_gnm
        self.current_f = self.open_local_or_remote_gnm(current_gnm,first_seek,size)
        #super(VirtualFile, self).__init__(self.current_f)

    def read(self,size):
        import sys
        sys.stderr.write('read called with ' + self.current_gnm + ' with size = ' + str(size) + '\n')
        #current_f = self.open_local_or_remote_gnm(self.current_gnm,self.first_seek)
        b = self.current_f.read(size)
        while len(b) < size:
            self.current_f.close()
            if self.current_gnm in genome_succession:
                self.current_gnm = genome_succession[self.current_gnm]
                self.current_f = self.open_local_or_remote_gnm(self.current_gnm,0,size)
                b += self.current_f.read(size-len(b))
            else:
                break
        return b

    def open_local_or_remote_gnm(self,gnm,seek_to,size):
        if protocol_matcher.match(gnm):
            return urllib2.urlopen(urllib2.Request(gnm, headers = {'Range': 'bytes='+str(seek_to)+'-'+str(seek_to+size-1)}))
        else:
            import io
            f = io.open(gnm,'rb')
            import sys
            sys.stderr.write('seeking to ' + str(seek_to) + '\n')
            f.seek(seek_to)
            return f

@app.route('/virtualgenome.fa.fai')
@app.route('/virtualgenome.fa.gz.fai')
def return_fai():
    return get_fai()

@app.route('/virtualgenome.fa.gz.gzi')
def return_gzi():
    return get_gzi()

@app.route('/virtualgenome.fa', methods=['HEAD'])
@app.route('/virtualgenome.fa.gz', methods=['HEAD'])
def satisfy_head_check():
    #FIXME: probably ought to include Response headers
    r = make_response()
    #r['Accept'] = 'Bytes';
    return r

@app.route('/virtualgenome.fa', methods=['GET'])
@app.route('/virtualgenome.fa.gz', methods=['GET'])
def return_fa_slice():
    import sys
    #FIXME: do this with logging, though probably not really even needed anymore
    for header in request.headers.keys():
        sys.stderr.write(header + ' ' + request.headers.get(header) + '\n')

    byterange = request.headers.get('Range')
    if byterange:
        #FIXME: do this with logging, though probably not really even needed anymore
        sys.stderr.write('byterange='+str(byterange)+'\n')
        return get_fa_slice(str(byterange))

def get_fai():
    """constuct a virtual genome by aggregating all the individual genomic fais
    if there is none yet for the current application context.
    """
    global virtual_fai
    if not virtual_fai:
        virtual_fai = create_virtual_fai()
    #FIXME: eventually, probably just write this to disk when constructed and return that way
    import sys
    sys.stderr.write(virtual_fai)
    return virtual_fai

def create_virtual_fai():
    """implementation details of virtual fai construction
    """
    global genome_sizes
    global genome_succession
    retval = ''
    offset = 0
    global virtual_fai_offsets
    last_genome = None;
    with open('config.txt') as cfg:
        genomes = cfg.read().splitlines()
    for gnm in (genomes):
        if gnm.startswith('#'):
            continue
        if protocol_matcher.match(gnm):
            thisfile = urllib2.urlopen(urllib2.Request(gnm+'.fai'))
            for line in thisfile: 
                data = line.split('\t')
                new_offset = int(data[2]) + offset
                data[2] = str(new_offset)
                retval += '\t'.join(data)
        else:
        #FIXME: DRY; should refactor logic that is same after streams opened
            with open(gnm+'.fai') as thisfile:
                for line in thisfile: 
                    data = line.split('\t')
                    new_offset = int(data[2]) + offset
                    data[2] = str(new_offset)
                    retval += '\t'.join(data)
        #figure out offset (ie what the concatenated offset would be in the virtual genome)
        lastseq_offset = int(data[2])
        lastseq_len = int(data[1])
        lastseq_line_len = int(data[3])
        lastseq_linedelim_len = int(data[4])-int(data[3])
        import math
        new_offset = offset + lastseq_offset + lastseq_len + lastseq_linedelim_len*(int(math.ceil(float(lastseq_len)/lastseq_line_len)))
        virtual_fai_offsets[offset:new_offset] = {'genome':gnm, 'offset':offset}
        offset = new_offset
        import sys
        sys.stderr.write('genome_size for ' + gnm + ' is ' + str(offset));
        genome_sizes[gnm] = offset
        #for VirtualFile strategy
        if last_genome is not None:
            genome_succession[last_genome] = gnm
        last_genome = gnm
    return retval

def get_gzi():
    """constuct a virtual genome by aggregating all the individual genomic gzis
    if there is none yet for the current application context.
    """
    global virtual_gzi
    if not genome_sizes:
        get_fai()
    if not virtual_gzi:
        virtual_gzi = create_virtual_gzi()
    return virtual_gzi

def create_virtual_gzi():
    """implementation details of virtual gzi construction
    """
    retval = ''
    start_compressed_offset = 0
    start_uncompressed_offset = 0
    total_entries = 0
    global virtual_gzi_offsets
    with open('config.txt') as cfg:
        genomes = cfg.read().splitlines()
    for gnm in genomes:
        if gnm.startswith('#'):
            continue
        #TODO: what if some genomes are not compressed? we probably can't mix and match
        if protocol_matcher.match(gnm):
            thisfile = urllib2.urlopen(urllib2.Request(gnm+'.gzi'))
            genome_size = genome_sizes[gnm]
            num_entries = unpack('<Q',thisfile.read(8))[0]
            total_entries += num_entries
            i = 0
            while i < num_entries:
                compressed_offset = unpack('<Q',thisfile.read(8))[0]
                uncompressed_offset = unpack('<Q',thisfile.read(8))[0]
                new_compressed_offset = compressed_offset + start_compressed_offset
                new_uncompressed_offset = uncompressed_offset + start_uncompressed_offset
                retval += pack('<Q',new_compressed_offset)
                retval += pack('<Q',new_uncompressed_offset)
                i+=1
        #FIXME: DRY; should refactor logic that is same after streams opened
        else:
            with open(gnm+'.gzi') as thisfile:
                genome_size = genome_sizes[gnm]
                num_entries = unpack('<Q',thisfile.read(8))[0]
                total_entries += num_entries
                i = 0
                while i < num_entries:
                    compressed_offset = unpack('<Q',thisfile.read(8))[0]
                    uncompressed_offset = unpack('<Q',thisfile.read(8))[0]
                    new_compressed_offset = compressed_offset + start_compressed_offset
                    new_uncompressed_offset = uncompressed_offset + start_uncompressed_offset
                    retval += pack('<Q',new_compressed_offset)
                    retval += pack('<Q',new_uncompressed_offset)
                    i+=1
        new_uncompressed_offset = start_uncompressed_offset + genome_sizes[gnm]
        #not sure using genome_size makes sense for compressed, but can't think of anything better
        #FIXME: follow IGV logic in determining end range which is some function of a fixed buffer size?
        new_compressed_offset = start_compressed_offset + genome_sizes[gnm]
        import sys
        sys.stderr.write('offset for ' + gnm + ' is ' + str(start_compressed_offset) + '\n')
        #the -1 was added to this after a commit that took place before a long hiatus from this code; not sure if it's a good thing but leaving for now while I try the VirtualFile strategy
        virtual_gzi_offsets[start_compressed_offset:new_compressed_offset-1] = {'genome':gnm, 'offset':start_compressed_offset}
        start_uncompressed_offset = new_uncompressed_offset
        start_compressed_offset = new_compressed_offset
    return pack('<Q',total_entries)+retval

def get_fa_slice(byterange):
    m = byterange_matcher.match(byterange)
    start = int(m.group(1))
    stop = int(m.group(2))
    if len(virtual_gzi_offsets):
        import sys
        sys.stderr.write('interval overlaps='+str(len(virtual_gzi_offsets[start]))+'\n')
        d = virtual_gzi_offsets[start].pop().data
    else:
        d = virtual_fai_offsets[start:stop].pop().data
    gnm = d['genome']
    import sys
    sys.stderr.write('byte range from ' + str(start) + ' to ' + str(stop) + ' was matched to genome ' + gnm + ' with offset ' + str(d['offset']) + '\n')
    f = VirtualFile(gnm,start-d['offset'],8192)
    r = make_response(send_file(f, 'X-application-fasta'), 206)
    r.headers['Content-Length'] = stop-start+1
    r.headers['Accept-Ranges'] = 'bytes'
    return r

