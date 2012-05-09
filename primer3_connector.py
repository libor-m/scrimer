
# A connector to the primer3_core executable.
# (it is not possible to keep a persistent primer3 process
#  using subprocess module - communicate() terminates the input
#  stream and waits for the process to finish)
#
# Author: Libor Morkovsky 2012
#
#TODO: process_records() taking a list of dicts, returning list of
# dicts using a single call of primer3 (should be more effective
# than to call primer3 for each record)
class BoulderIO:
    @classmethod
    def parse(self, string):
        """parse a BoulderIO string (KEY=VAL\n)
        return a list of records, where each record is a dictionary
        end of the string implies a single '=\n' (record separator)
        """
        record_strings = string.split("=\n")
        return [dict(tuple(line.split("=", 1)) for line in record.split("\n") if len(line) > 3) for record in record_strings if len(record) > 3]

    @classmethod
    def deparse(self, records):
        """accept a dict or a list of dicts
        produce a BoulderIO string (KEY=VAL\n)
        with records separated by '=\n' when the input was a list
        """
        
        # unify the input, create a list with single element
        if type(records) == dict:
            records = [records]

        return "\n=\n".join("\n".join("=".join(kval) for kval in record.iteritems()) for record in records) + "\n=\n"


class Primer3:
    def __init__(self, p3path="primer3_core", **kwargs):
        # store path to primer3
        self.p3path = p3path
        
        # add stringized versions of all kwargs to default args
        self.default_params = {}
        str_kw = dict((key, str(val)) for key, val in kwargs.iteritems())
        self.default_params.update(str_kw)

    def process_record(self, record):
        """OBSOLETE, superseded by call()
        Merge the record with default_params, the record taking precedence,
        call the primer3 worker, 
        parse the output and return a list of dictionaries 
        {RIGHT:[], LEFT:[], PAIR:[], INTERNAL:[]} for each input record
        """
        
        # merge the defaults with current query
        full_record = dict(self.default_params.items() + record.items())
        
        # call primer3
        import subprocess
        self.child = subprocess.Popen([self.p3path, '-strict_tags'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = self.child.communicate(BoulderIO.deparse(full_record))
        
        #TODO: check for errors in stderr
        
        #TODO: process multiple BoulderIO records in the output
        results = BoulderIO.parse(out)[0]
        
        # parse the results to {RIGHT:[], LEFT:[], PAIR:[], INTERNAL:[]}
        sides = ['RIGHT', 'LEFT', 'PAIR', 'INTERNAL']
        # add a separate list for each side, .fromkeys() adds a reference to single instance
        primers = dict((side, []) for side in sides)
        for side in sides:
            nret_key = 'PRIMER_%s_NUM_RETURNED' % side
            nret = int(results.get(nret_key, 0))
            
            # extract the values for each single primer and put those to 
            # equivalent key
            for num in xrange(nret):
                template = 'PRIMER_%s_%d_' % (side, num)
                primer_keys = filter(lambda k: template in k, results.iterkeys())
                primer = dict((key[len(template):], results[key]) for key in primer_keys)
                
                primers[side].append(primer)

        return primers

    def call(self, records):
        """Merge each of the records with default_params, the record taking precedence,
        call the primer3 worker, 
        parse the output and return a list of dictionaries,
        {RIGHT:[], LEFT:[], PAIR:[], INTERNAL:[]} for each input record
        uppercase keys (in the result) are the original names from BoulderIO format,
        lowercase keys have no direct equivalent (position, other-keys)
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
