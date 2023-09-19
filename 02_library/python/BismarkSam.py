class BismarkSam(object):
    """
    Class to import and manipulate the XM tags of SAM files exported from Bismark.
    """
    def __init__(self, read) -> None:
        """
        Initialise class function
        """
        if read[13][:5] != "XM:Z:":
            raise ValueError("Element 13 of this read is not an XM tag. This read does not appear to be a Bismark-sorted SAM file.")
        self.id = read[0]
        self.chr = read[2]
        self.length = len(read[9])
        self.xm_tag = read[13][5:]#self.trim_XM_tag(read[13][5:])
    
    def trim_XM_tag(self, xm_tag):
        """
        Pull out the XM tag from a read in a SAM file and trim away the non-cytosines
        """
        trim_tag = [c for c in xm_tag if (c.islower() or c.isupper()) ]
        
        return trim_tag
    
    def count_mC(self):
        """
        Count how many cytosines on a read are methylated, unmethylated, and the total.
        """
        upper = 0
        lower = 0
        up="HXZ"
        lo="hxz"
        for i in self.xm_tag:
            if i in up:
                upper+=1
            elif i in lo:
                lower+=1
        return [upper, lower]

    def mC_cluster(self) :
        """
        Check whether all the methylated cytosines appear together in a read
        (there  are no unmethylated cytosines between methylated cytosines).
        
        It is important that non-cytosines have been removed from the read.
        """
        flag = False
        index = 0
        n = len(self.xm_tag)
        
        # Check for clusters in reads with two or more cytosines.
        while index < n:
            if self.xm_tag[index].isupper():
                if (flag == True) :
                    return False
                while index < n and self.xm_tag[index].isupper():
                    index += 1
                flag = True
            else :
                index += 1
        return True
    
    def mC_per_read(self):
        """
        Summarise the number of methylated reads, unmethylated reads, total read
        length, and whether or not cytosines are occur next to one another
        """
        total = self.length
                        
        mC, uC = self.count_mC()
        if mC > 1:
            cluster = self.mC_cluster()
        else:
            cluster = 'NA'
        
        return [self.chr, mC, uC, total, cluster]