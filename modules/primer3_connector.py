"""
A Python interface to the primer3_core executable.

    TODO: it is not possible to keep a persistent primer3 process
     using subprocess module - communicate() terminates the input
     stream and waits for the process to finish

Author: Libor Morkovsky 2012
"""

# This file is a part of Scrimer.
# See LICENSE.txt for details on licensing.
#    Copyright (C) 2012, 2013 Libor Morkovsky

class BoulderIO:
    """Provides Python interface for ``BoulderIO`` format used by Primer3.
    """

    @classmethod
    def parse(self, string):
        r"""Parse a BoulderIO string ``(KEY=VAL\n)``
        return a list of records, where each record is a dictionary
        end of the string implies a single ``'=\n'`` (record separator).
        """
        record_strings = string.split("=\n")
        return [dict(tuple(line.split("=", 1)) for line in record.split("\n") if len(line) > 3) for record in record_strings if len(record) > 3]

    @classmethod
    def deparse(self, records):
        r"""Accepts a dict or a list of dicts, produces a BoulderIO string ``(KEY=VAL\n)``
        with records separated by ``'=\n'``.
        """
        
        # unify the input, create a list with single element
        if type(records) == dict:
            records = [records]

        return "\n=\n".join("\n".join("=".join(kval) for kval in record.iteritems()) for record in records) + "\n=\n"


class Primer3:
    """Wraps Primer3 executable. `kwargs` are converted to strings and used as default parameters
    for each call of primer3 binary.
    """
    def __init__(self, p3path="primer3_core", **kwargs):
        # store path to primer3
        self.p3path = p3path
        
        # add stringized versions of all kwargs to default args
        self.default_params = {}
        str_kw = dict((key, str(val)) for key, val in kwargs.iteritems())
        self.default_params.update(str_kw)

    def call(self, records):
        """Merge each of the records with `default_params`, the record taking precedence,
        call the ``primer3`` binary, 
        parse the output and return a list of dictionaries,
        ``{RIGHT:[], LEFT:[], PAIR:[], INTERNAL:[]}`` for each input record
        uppercase keys (in the result) are the original names from BoulderIO format,
        lowercase keys have no direct equivalent in primer3 output (``position``, ``other-keys``)
        """

        # merge the defaults with current query
        full_records = [dict(self.default_params.items() + record.items()) for record in records]

        # call primer3
        import subprocess
        self.child = subprocess.Popen([self.p3path, '-strict_tags'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = self.child.communicate(BoulderIO.deparse(full_records))
        
        # simple check for errors in stderr
        if len(err):
            raise Exception(err)
        
        results = BoulderIO.parse(out)
        
        # parse the results to {RIGHT:[], LEFT:[], PAIR:[], INTERNAL:[]}
        sides = ['RIGHT', 'LEFT', 'PAIR', 'INTERNAL']
        primers = []
        
        for result in results:
            # primers for current result
            res_primers = dict((side, []) for side in sides)
            used_keys = []
            
            for side in sides:
                nret_key = 'PRIMER_%s_NUM_RETURNED' % side
                nret = int(result.get(nret_key, 0))
                used_keys.append(nret_key)
                
                # extract the values for each single primer and put those to 
                # equivalent key
                for num in xrange(nret):
                    template = 'PRIMER_%s_%d_' % (side, num)
                    primer_keys = filter(lambda k: template in k, result.iterkeys())
                    primer = dict((key[len(template):], result[key]) for key in primer_keys)
                    
                    # extract the position, which itself has no extractible name in BoulderIO
                    # only 'PRIMER_LEFT_0'
                    if side != 'PAIR':
                        pos_key = template[:len(template)-1]
                        primer['position'] = result.get(pos_key, "#error!")
                        used_keys.append(pos_key)
                    
                    # keep track of keys used in current record
                    used_keys.extend(primer_keys)
                    
                    res_primers[side].append(primer)
            
            # store all the unused keys for current result
            res_primers['other-keys'] = dict((key, result[key]) for key in result.iterkeys() if key not in used_keys)
            primers.append(res_primers)
        
        return primers

if __name__ == "__main__":
    print "Running tests"
    
    import textwrap
    
    record = BoulderIO.parse(textwrap.dedent(
    """
    SEQUENCE_ID=example
    SEQUENCE_TEMPLATE=GTAGTCAGTAGACGATGACTACTGACGATGCAGACNACACACACACACACAGCACACAGGTATTAGTGGGCCATTCGATCCCGACCCAAATCGATAGCTACGATGACG
    SEQUENCE_TARGET=37,21
    PRIMER_PICK_INTERNAL_OLIGO=0
    PRIMER_OPT_SIZE=18
    PRIMER_MIN_SIZE=15
    PRIMER_MAX_SIZE=21
    PRIMER_MAX_NS_ACCEPTED=3
    PRIMER_PRODUCT_SIZE_RANGE=50-100
    """))

    record_no_res = BoulderIO.parse(textwrap.dedent(
    """
    SEQUENCE_ID=example
    SEQUENCE_TEMPLATE=GTAGTCAGTAGACNATGACNACTGACGATGCAGACNACACACACACACACAGCACACAGGTATTAGTGGGCCATTCGATCCCGACCCAAATCGATAGCTACGATGACG
    SEQUENCE_TARGET=37,21
    PRIMER_TASK=pick_detection_primers
    PRIMER_PICK_LEFT_PRIMER=1
    PRIMER_PICK_INTERNAL_OLIGO=1
    PRIMER_PICK_RIGHT_PRIMER=1
    PRIMER_OPT_SIZE=18
    PRIMER_MIN_SIZE=15
    PRIMER_MAX_SIZE=21
    PRIMER_MAX_NS_ACCEPTED=1
    PRIMER_PRODUCT_SIZE_RANGE=75-100
    SEQUENCE_INTERNAL_EXCLUDED_REGION=37,21
    """))

    default_params = BoulderIO.parse(textwrap.dedent(
    """
    PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/opt/primer3/bin/primer3_config/
    PRIMER_MAX_NS_ACCEPTED=0
    PRIMER_EXPLAIN_FLAG=1
    """))[0]

    print "Testing BoulderIO, single record:", 
    record_dp = BoulderIO.deparse(record)
    record_reparsed = BoulderIO.parse(record_dp)
    if record == record_reparsed:
        print "OK"
    else:
        print "Failed!"

    print "Testing BoulderIO, two records:", 
    two_records = record + record_no_res
    
    record_dp = BoulderIO.deparse(two_records)
    record_reparsed = BoulderIO.parse(record_dp)
    if two_records == record_reparsed:
        print "OK"
    else:
        print "Failed!"
    
    print "Testing Primer3, single record:", 
    p3 = Primer3(**default_params)

    # test for single record
    res = p3.call(record)
    if res[0]['RIGHT'][0]['SEQUENCE'] == 'GTCGGGATCGAATGGCCC':
        print "OK"
    else:
        print "Failed!"
 
    # test for multiple records
    print "Testing Primer3, two records:", 
    res = p3.call(two_records)

    # second record should produce no results
    if len(res[1]['RIGHT']) == 0:
        print "OK"
    else:
        print "Failed!"
    
    # if no exception occurs, the test should be OK
    print "Tests ran OK"
