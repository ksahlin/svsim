import os
import subprocess
import tempfile
from numpy.random import lognormal
from pkg_resources import resource_filename

from svsim.util import calculate_num_reads, get_genome_length
from svsim.reads.isim import IReadSimulator

##
# Simulates reads using metasim.
#
class LogNSimulator( IReadSimulator ):
    ##
    # Constructor.
    #
    def __init__(self):
        if self.mean < 10 and  1 < self.std <  2: # Well defined distribution
            pass
        else: # user probably forgot that mu and std needs to be specified in log base
            raise ValueError( "mu and std needs to be specified in log base (usually mu < 10 and 1 < sigma < 2)" )


    ##
    # @see IReadSimulator.simulate
    #
    def simulate(self, genome_path, output_prefix):
        if not self.read_length == 100:
            raise ValueError( "Read length must be 100 for metasim." )

        genome_sequence = get_genome_sequence( genome_path )
        genome_length = len(genome_sequence)
        dna_library = DNAseq(self.read_length, self.coverage, mean=self.mean, stddev=self.std, distribution='lognormal')
        dna_library.simulate_pe_reads(genome_sequence)

        reads1 = open(os.path.join(output_prefix, "_pe1.fa"),'w')
        reads2 = open(os.path.join(output_prefix, "_pe2.fa"),'w')
        i=0
        for read in dna_library.fasta_format():
            if i%2==0:
                reads1.write(read)
            else:
                reads2.write(read)
            i+=1

        # output_dir = tempfile.mkdtemp( )

        # subprocess.call( [ "MetaSim", "cmd",
        #                    "-r", str( num_reads ), 
        #                    "-m",
        #                    "-g", self.error_model,
        #                    "-2", self.error_model,
        #                    "--empirical-pe-probability", "100",
        #                    "--clones-mean", str( self.mean ),
        #                    "--clones-param2", str( self.std ),
        #                    "-d", output_dir,
        #                    genome_path ], stdout = open( "/dev/null", "w" ) )


class PairedEndRead(object):
    """docstring for PairedEndRead"""
    def __init__(self,distribution = 'normal',mean=None,sigma=None,read_length= None, min_size=None, max_size=None):
        super(PairedEndRead, self).__init__()
        self.distribution = distribution
        self.mean = mean
        self.sigma = sigma
        self.read_length = read_length
        self.min_size = min_size
        self.max_size = max_size

    def generate(self, reference_sequence, read_index):
        if self.distribution == 'normal':
            self.fragment_length = int(random.gauss(self.mean,self.sigma))
        elif self.distribution == 'uniform':
            self.fragment_length = int(random.uniform(self.min_size,self.max_size))
        elif self.distribution == 'lognormal':
            self.fragment_length = int(lognormal(self.mean, self.sigma)[0]) # one sample at a time to conform with the implementation...

        if self.fragment_length >= len(reference_sequence): 
            raise Exception("To short reference sequence length for \
                simulated read. \nRead fragment: {0}\nTranscript \
                length:{1}".format(self.fragment_length,len(reference_sequence)))
        
        self.start_pos = random.randrange(len(reference_sequence) - self.fragment_length)
        self.read1 = reference_sequence[self.start_pos : self.start_pos + self.read_length]
        self.read2 = reverse_complement(reference_sequence[self.start_pos + self.fragment_length - self.read_length : self.start_pos+self.fragment_length])
        self.reference_accession = "reference_genome"
        self.read_index = read_index

    def fastq_format(self):
        r1= '@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{2}\
        \n@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{3}\n'.format(self.read_index,
            self.reference_accession,self.read1,'J'*self.read_length)
        r2= '@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{2}\
        \n@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{3}\n'.format(self.read_index,
            self.reference_accession,self.read2,'J'*self.read_length)        
        yield r1
        yield r2

    def fasta_format(self):
        r1= '>HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{2}\n'.format(self.read_index,
            self.reference_accession,self.read1)
        r2= '>HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{2}\n'.format(self.read_index,
            self.reference_accession,self.read2)        
        yield r1
        yield r2


class DNAseq(object):
    """docstring for DNAseq"""
    def __init__(self,read_length, coverage, mean=None,stddev=None, min_size=None, max_size = None, distribution='normal'):
        super(DNAseq, self).__init__()
        self.distribution = distribution
        
        if self.distribution == 'normal' or self.distribution == "lognormal":
            self.mean = mean
            self.stddev = stddev
        elif self.distribution == 'uniform':
            self.min_size = min_size
            self.max_size = max_size

        self.read_length = read_length
        self.coverage = coverage

    def simulate_pe_reads(self, genome_sequence):
        """
        Arguments:
        """
        genome_length = len(genome_sequence)
        number_of_reads = calculate_num_reads(self.coverage, self.read_length, genome_length)  # Specifies the number of simulated read pairs (related to insertion size length of genome and coverage
    
        self.reads = []
        
        
        for i in range(number_of_reads):
            if self.distribution == 'normal':
                read_pair = PairedEndRead(distribution=self.distribution, mean=self.mean, sigma=self.stddev, read_length=self.read_length)
            elif self.distribution == 'lognormal':
                read_pair = PairedEndRead(distribution=self.distribution, mean=self.mean, sigma=self.stddev, read_length=self.read_length)
            elif self.distribution == 'uniform':
                read_pair = PairedEndRead(distribution=self.distribution, min_size=self.min_size,max_size=self.max_size,read_length=self.read_length)
            
            read_pair.generate(genome_sequence, i)
            self.reads.append(read_pair)


    def fastq_format(self):
        for pe_read in self.reads:
            for mate in pe_read.fastq_format():
                yield mate

    def fasta_format(self):
        for pe_read in self.reads:
            for mate in pe_read.fasta_format():
                yield mate