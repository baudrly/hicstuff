# -*- coding: utf-8 -*-

import numpy as np
import warnings
import hicstuff as hcs
import io
from collections import OrderedDict,Counter
import types

try:
    basestring
except NameError:
    basestring = str

class HiCMatrix(object):
    """Constructs an Hi-C matrix object from various data formats.
    Annotations, genomic information (positions, contigs...) are automatically
    added as needed.
    """
    
    def __init__(self, matrix_file, name="", size="auto", n_contigs=1,
                 genome=[""], annotations={}, positions="auto", contigs="auto", 
                level=0, GC="auto", sparsity="auto", structure=[]):
                    
        """Loads the Hi-C instance into attributes as follows:
        -matrix_file: matrix-like object such as an array-like or a path to a 
        matrix file (in csv or Dade format).
        E.g.: '/home/user/matrix.txt', np.array([0]), [[1,2],[2,5]], etc.
        
        -name: name of the species or any kind of suitable label.
        E.g.: 'Saccaromyces cerevisiae', 'test_3', etc.
        
        -size: matrix size (in pixels, or fragments, or bins). Matrix length
        by default.
        
        -n_contigs: number of contigs or chromosomes. Defaults to 1.
        
        -genome: a list (of same length as n_contigs) of Biopython Seq records 
        or strings, or a generator iterating on either, or a file object 
        containing the sequence, or a file path to the genome.
        E.g.: '~/genomes/test.fasta', '/home/user/test2.gbk', ['AATTC','CWN'],
        [<class 'Bio.SeqRecord.SeqRecord'>], etc.
        
        -annotations: a dictionary in the format {'annotation_name':annotation}
        where annotation is an array-like of same length as the size attribute.
        E.g.: {'annotation1':[0,0,0,5,0],'annotation2':[1,2,-1.5,1,1]}, {}, etc.
        
        -positions: an array-like (of same length as size) containing the 
        starting positions (in bp) of each pixel (or fragment, or bin) in the 
        matrix within their respective contig.
        E.g.: [[0,10000,20000,30000,35864],[0,10000,20000]], 
        [[0,8756,15422,28550,35788]], etc.
        If positions is set to "auto", it will be inferred from matrix 
        construction (assuming the input contains position information) before
        being set to default.
        If input is an integer, positions will be generated using that integer's value
        as binning across the genome.
        
        -contigs: an array-like (of same length as size) containing integers
        indicating what contig each pixel (or fragment, or bin) in the matrix
        belongs to. All ones by default.
        E.g.: [1,1,1,1,1,1], [1,1,1,2,2,3,3,3], etc.
        If contigs is set to "auto", it will be inferred from matrix 
        construction (assuming the input contains position information) before
        being set to default.
        
        -level: an integer indicating how the matrix is binned. 
        If level < 50, it means the matrix is fragment-binned, i.e. each pixel
        in an n-level matrix is equal to the sum of the pixels in the (n-1)-level
        corresponding matrix contained in a square surrounding that pixel, e.g.
        if level 0 is [[1,2,3,4,5,6],
                       [7,8,9,0,1,2],
                       [3,4,5,6,7,8],
                       [9,0,1,2,3,4],
                       [5,6,7,8,9,0],
                       [1,2,3,4,5,6]]
        level 1 will be:
                       [[38,39],
                        [34,41]]
        In other words, fragments are binned into groups of 3 such that each bin
        from level n represents 3^n fragments from the initial fragment list.
        If level is 50 or more, it means the matrix is kb-binned, i.e. signal
        is summed and distributed along bins of equal length (equal to level,
        so greater than 50 pb).
        
        -GC: an array-like (of same length as size) containing, 
        for each contig, the GC content of each bin along the genome in 
        percentage form. 
        E.g.: [0,50,40,60], etc.
        
        -Sparsity: an array-like (of same length as size_contig) containing 
        contact sparsity (i.e. sum of contacts) for each pixel (or fragment,
        or bin).
        """
        #class attribute
        KB_BINNING = 10000.0
        
        #First load any specifically indicated attributes
        self.name = name
        self.positions = positions
        self.contigs = contigs
        self.n_contigs = n_contigs
        self.level = level
        self.structure = structure
        #Initialize matrix and relevant coordinate information
        self._init_matrix(matrix_file)
        self.size = size if not size == "auto" else len(self.matrix)

        #Then initialize genome and annotations
        self.genome = genome
        self._init_genome(genome)
        self.annotations = annotations
        self._init_annotations(annotations)
       
        #If positions and contigs weren't specified at first or inferred
        #from previous initialization, set some sane default values
        if self.positions == "auto":
            self.positions = [int(j) for i in range(self.n_contigs) for j in np.linspace(0, self.size*KB_BINNING, self.size)]
        elif isinstance(self.positions, (int,float,np.int32,np.int64,np.float32,np.float64)):
            binning = self.positions
            self.positions = [int(j) for i in range(self.n_contigs) for j in np.linspace(0, self.size*binning, self.size)]
            self.level = binning
            
        if self.contigs == "auto":
            self.contigs = np.ones(self.size)
        
        #GC and sparsity are used all the time so they deserve attributes on
        #their own
        self.GC = GC
        if GC == "auto" and genome[-1]:
            self.load_GC(genome)
            
        self.sparsity = sparsity
        if sparsity == "auto":
            self.load_sparsity()
    
    def _init_matrix(self, matrix_file):
        #Check for various data types before initializng matrix
    
        #Path:
        matrix = matrix_file
        
        if isinstance(matrix_file, basestring):
            try:
                matrix = np.genfromtxt(matrix_file, dtype=None, filling_values=0)
                print("test")
            except ValueError as e:
                print("Input must be matrix or path to matrix.")
                print(str(e))
                print("test3")
        #Any kind of array that isn't a numpy array:
        elif not type(matrix_file) is np.ndarray:
            matrix = np.array(matrix_file)
        if len(matrix.shape) < 2 and len(matrix) > 1:
            raise ValueError("Input matrix must be sparse or 2D")
        elif len(matrix.shape) == len(matrix) == 1:
            matrix.shape = (1,1)
        is_sparse = min(matrix.shape) == 3 and max(matrix.shape) != 3
        is_horizontal = is_sparse and matrix.shape[0] > matrix.shape[1]
        is_square = not is_sparse and matrix.shape[0] == matrix.shape[1]
        if not is_square:
            warnings.warn("Input matrix was found not to be square. I will try to sort out the extra lines, but things may be messy.", RuntimeWarning)
        
        if is_horizontal:
            matrix = matrix.T
        
        #Check matrix type and intialize accordingly
        if matrix[0,:].dtype.type is np.str_ or matrix[0,:].dtype.type is np.bytes_:
            if matrix[:,0].dtype.type is np.str_ or matrix[:,0].dtype.type is np.bytes_:
                self._init_from_dade_matrix(matrix)               
            else:
                self._init_from_header_matrix(matrix)
        else:
            self.matrix = matrix  
        
    def _init_from_dade_matrix(self, matrix, *args, **kwargs):
        """Extracts matrix, positions and contigs from input Dade matrix."""
        self.matrix, header = np.array(matrix[1:,1:],dtype=int), matrix[0,:]
        parsed_header = list(zip(*[str(h)[:-1].strip('"').strip("'").split("~") for h in header[1:]]))
        if len(parsed_header) == 2:
            contigs, positions =  parsed_header
        else:
            global_index, contigs, index, positions, _ = parsed_header
        
        #yields somethings like {'subtilis_chr1': 1, 'vibrio_plasmid1:'2'}, and so on
        contig_map = dict(zip(OrderedDict.fromkeys(contigs),np.array(range(1,len(set(contigs))+1))))
        
        positions = np.array(positions,dtype=int)
        indices_zero = np.arange(len(positions))[positions==0]
        self.contigs = np.array([contig_map[contig] for contig in contigs])
        self.positions = np.split(positions, indices_zero)
        
    def _init_from_header_matrix(self, matrix, *args, **kwargs):
        """Removes headers from matrix file and loads the resulting matrix
        into the matrix attribute."""
        self.matrix, _ = matrix[1:,], matrix[0,:]
        
    def _init_from_sparse_coo_matrix(self, sparse, *args, **kwargs):
        
        if sparse[0].dtype is np.str_:
            sparse = sparse[1:,:]
        try: 
            from scipy.sparse import coo_matrix
            matrix = coo_matrix(sparse[:,2],(sparse[:,0],sparse[:,1]))
        except ImportError as e:
            print("I can't load a sparse matrix.")
            print(str(e))
        self.matrix = matrix
        
    def _init_GC(self):
        self.load_GC(self.genome)
        
    def load_GC(self,genome):
        """Compute GC information based on currently loaded genome and
        store it in GC attribute."""
        GC = []
        for contig in genome:
            for i in range(len(self.positions)):
                portion = contig[self.positions[i]:self.positions[i+1]]
                try:
                    from Bio import SeqUtils
                    gc = SeqUtils.GC(portion)
                except ImportError: #Fallback if Biopython isn't available
                    gc = self.partial_GC(portion)
                GC.append(gc)
                
        self.GC = np.array(GC)
        
    def load_sparsity(self):
        self.sparsity = self.matrix.sum(axis=0)
        
    def partial_GC(self,portion):
        """Manually compute GC content percentage in a DNA string, taking
        ambiguous values into account (according to standard IUPAC notation)."""
        sequence_count = Counter(portion)
        gc = sum([sequence_count[i] for i in 'gGcCsS'])
        + sum([sequence_count[i] for i in 'DdHh'])/3.0
        + 2*sum([sequence_count[i] for i in 'VvBb'])/3.0
        + sum([sequence_count[i] for i in 'NnYyRrKkMm'])/2.0
        return 0 or 100*gc
                
    def _init_genome(self, genome, load_full=True):
        self.load_genome(genome, load_full)
    
    def load_genome(self, genome, load_full=True):
        """Load a genome into the HiC object. Requires Biopython.
        Genome can take different formats:
        -Path to genome file, genbank or fasta (by default)
        -File object
        -Genome generator (typical output of Biopython parsing)
        
        If load_full is False, the genome attribute will contain a generator
        pointing to the genome's contigs. Otherwise, the entirety of the genome
        will be loaded into the genome attribute.
        """
        
        try:
            from Bio import SeqIO
            if isinstance(genome, basestring):
                print("I detected a genome file path.")
                if genome.split(".")[-1] in ["gbk","genbank"]:
                    genome_format = "genbank"
                else:
                    genome_format = "fasta"
                genome = open(genome,'r')
            if isinstance(genome, io.IOBase):
                print("I detected a genome file.")
                handle = genome
                genome = SeqIO.parse(handle,genome_format)
                
            if load_full and isinstance(genome, types.GeneratorType):
                self.genome = list(genome)
                if len(self.genome) != len(set(self.contigs)):
                    warnings.warn("Warning, genome and contig information mismatch.", RuntimeWarning)
            else:
                self.genome = genome
                
        except ImportError as e:
            print("I couldn't import Biopython which is needed to handle genomes.")
            print(str(e))
        except IOError as e:
            print("I couldn't load input genome file.")
            print(str(e))
    
    def _init_annotations(self, annotations):
        for name,annotation in annotations.items():
            self.load_annotation(annotation, name)
                
    def load_annotation(self, annotation_file, name):
        """Store any annotation (or file path pointing to one) in the form of
        a dictionary {name: annotation} in the annotation attribute.
        """
        annotation = annotation_file
        if isinstance(annotation, basestring):
            annotation = np.genfromtxt(annotation_file,dtype=None)           
        self.annotations[name] = annotation
        
    def get_contig_borders(self):
        borders = [0]
        n = len(self.contigs)
        for i in range(1,n):
            if self.contigs[i-1] != self.contigs[i]:
                borders.append(i)
        borders.append(n-1)
        return np.array(borders)
        
        
    def __getitem__(self,index):
        """If index is a slice, return a diagonal-centered square matrix of
        size equal to the slice's with all genome information and annotations
        being sliced accordingly.
        If index is an integer i, return the matrix (and relevant genome
        information) corresponding to the i-th contig.
        """
        if isinstance(index, slice):
            return self._get_slice(index)
        elif isinstance(index, (int,np.int64,np.int32)):
            return self._get_contig(index)
        elif isinstance(index, (tuple,list)):
            if len(index) == 1:
                return self._get_contig(index[0])
            else:
                new_matrix = self._get_contig(index[0])
                for i in index[1:]:
                    new_matrix = new_matrix + self._get_contig(i)
            return new_matrix
            
    def _get_contig(self,index):
        contig_borders = self.get_contig_borders()
        start,end = contig_borders[index],contig_borders[index+1]
        new_matrix = self.matrix[start:end,start:end]
        new_size = len(new_matrix)
        new_positions = self.positions[start:end]
        new_annotations = {name:annotation[start:end] for name,annotation in self.annotations.items()}
        new_contigs = self.contigs[start:end]
        new_n_contigs = len(set(new_contigs))
        try:
            new_sparsity = self.sparsity[start:end]
        except IndexError:
            new_sparsity = []
        try:
            new_GC = self.GC[start:end]
        except IndexError:
            new_GC = []
        try:
            new_genome = [self.genome[index]]
        except IndexError:
            new_genome = ['']
        new_structure = self.structure
        return HiCMatrix(new_matrix, self.name, new_size, new_n_contigs,
                         new_genome, new_annotations, new_positions, new_contigs,
                         self.level, new_GC, new_sparsity, new_structure)
                         
    def _get_slice(self,index):
        if index.stop >= index.start:
            new_matrix = self.matrix[index.start:index.stop:index.step,index.start:index.stop:index.step]
        else:
            new_matrix = np.array([0])
        new_size = len(new_matrix)
        new_positions = self.positions[index.start:index.stop:index.step]
        new_annotations = {name:annotation[index.start:index.stop:index.step] for name,annotation in self.annotations.items()}
        new_contigs = self.contigs[index.start:index.stop:index.step]
        new_n_contigs = len(set(self.contigs))
        new_sparsity = self.sparsity[index.start:index.stop:index.step]
        new_GC = self.GC[index.start:index.stop:index.step]
        
        #Since genome is an array of contigs, first take the (presumably)
        #incomplete end of the first contig, the (presumably) incomplete start
        #of the last one, and all the contigs inbetween.
        first_contig = int(self.contigs[index.start])
        last_contig = int(self.contigs[index.stop])
        if first_contig > last_contig:
            raise ValueError("Contigs are not well-ordered.")
        
        if isinstance(self.genome, list) and self.genome[0]:
            new_genome = [self.genome[first_contig][self.positions[index.start]:]]
            new_genome += [self.genome[contig] for contig in range(first_contig+1,last_contig)]
            new_genome += [self.genome[last_contig][:self.positions[index.stop]]]
        elif hasattr(self.genome, 'next'):
            new_genome = self.genome #no support for huge genomes for now
        else:
            new_genome = self.genome
                
        new_structure = self.structure #placeholder for now
        return HiCMatrix(new_matrix, self.name, new_size, new_n_contigs,
                         new_genome, new_annotations, new_positions, new_contigs,
                         self.level, new_GC, new_sparsity, new_structure)
    
    def __add__(self, x):
        """Performs matrix concatenation along with all relevant annotations.
        """
        new_size = self.size + x.size
        new_matrix = np.zeros((new_size,new_size))
        new_matrix[:self.size,:self.size] = self.matrix[:,:]
        new_matrix[self.size:,self.size:] = x.matrix[:,:]
        new_positions = np.array(list(self.positions)+list(x.positions))
        new_annotations = dict(self.annotations, **x.annotations)
        new_contigs =  np.array(list(self.contigs) + list(x.contigs))
        new_n_contigs = self.n_contigs + x.n_contigs
        new_sparsity = np.array(list(self.sparsity) + list(x.sparsity))
        new_GC = np.array(list(self.sparsity) + list(x.sparsity))
        new_genome = self.genome + x.genome
        new_structure = self.structure + x.structure
        return HiCMatrix(new_matrix, self.name, new_size, new_n_contigs,
                         new_genome, new_annotations, new_positions, new_contigs,
                         self.level, new_GC, new_sparsity, new_structure)
                         
    def load_structure(self,matrix,contigs):
        self.structure = hcs.tostruct(self.matrix)
        
    def write_pdb(self,filename=None):
        filename = self.name if filename is None else filename
        hcs.topdb(self.matrix,filename,self.contigs,self.annotations)
        print("Written to "+str(filename))
