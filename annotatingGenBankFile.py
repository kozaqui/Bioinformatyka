__author__ = 'DF from Bio Stars'
'''
https://www.biostars.org/p/57549/
Here's a small snippet of code for people (like myself) struggling through BioPython. I'm writing code to fix some annotations that are improperly annotated as being exact when they should really be open-ended (i.e. it's not certain where the feature actually starts). This bit was easy, but I got totally stuck for hours on actually adding the fixed annotation to a new or existing genbank file.

So here's a step-by-step guide (in runnable python code) to creating a new SeqRecord, annotating it with a SeqFeature, and then overwriting the SeqFeature in Biopython (I'm using python 2.7). Obviously this can be re-jigged to import & mess with existing genbank files. I wrote it for ease of understanding, not concise code, so sorry if it's a bit verbose. Simple for some people I know, but hopefully someone finds this useful.
'''



################ A: Make a SeqRecord ################

# 1. Create a sequence

from Bio.Seq import Seq
my_sequence = Seq("GATCGATCGATCGATCGATCGATCGATCGATC")

# 2. Create a SeqRecord and assign the sequence to it

from Bio.SeqRecord import SeqRecord
my_sequence_record = SeqRecord(my_sequence)

# 3. Assign an alphabet to the sequence (in this case DNA)

from Bio.Alphabet import generic_dna
my_sequence_record.seq.alphabet = generic_dna

# This is the minimum required info for BioPython to be able to output
# the SeqRecord in Genbank format.
# You probably would want to add other info (e.g. locus, organism, date etc)

#optional: print the SeqRecord to STDOUT in genbank format.. note there are no features on it yet.
print("\nThis bit is the SeqRecord, printed out in genbank format, with no features added.\n")
print(my_sequence_record.format("gb"))

################ B: Make a SeqFeature ################

# 1. Create a start location and end location for the feature
#    Obviously this can be AfterPosition, BeforePosition etc.,
#    to handle ambiguous or unknown positions

from Bio import SeqFeature
my_start_pos = SeqFeature.ExactPosition(2)
my_end_pos = SeqFeature.ExactPosition(6)

# 2. Use the locations do define a FeatureLocation
from Bio.SeqFeature import FeatureLocation
my_feature_location = FeatureLocation(my_start_pos,my_end_pos)

# 3. Define a feature type as a text string
#     (you can also just add the type when creating the SeqFeature)
my_feature_type = "CDS"

# 4. Create a SeqFeature
from Bio.SeqFeature import SeqFeature
my_feature = SeqFeature(my_feature_location,type=my_feature_type)

# 5. Append your newly created SeqFeature to your SeqRecord

my_sequence_record.features.append(my_feature)

#optional: print the SeqRecord to STDOUT in genbank format, with your new feature added.
print("\nThis bit is the SeqRecord, printed out in genbank format, with a feature added.\n")
print(my_sequence_record.format("gb"))

################ C: Overwrite an existing SeqFeature ################

# 1. Create a start location and end location for the feature..
#    This bit is obviously a repeat of "B: Make a SeqFeature" above,
#    normally I'd pull it out to a function, but I'm trying to be explicit here

from Bio import SeqFeature
my_start_pos = SeqFeature.ExactPosition(3)
my_end_pos = SeqFeature.ExactPosition(7)

# 2. Use the locations do define a FeatureLocation
from Bio.SeqFeature import FeatureLocation
my_feature_location2 = FeatureLocation(my_start_pos,my_end_pos)

# 3. Define a feature type as a text string
#    (or you can also just add the type when creating the SeqFeature)
my_feature_type2 = "ABC"

# 4. Create a SeqFeature
from Bio.SeqFeature import SeqFeature
my_feature2 = SeqFeature(my_feature_location2,type=my_feature_type2)

my_sequence_record.features[0]=my_feature2

#optional: print the SeqRecord to STDOUT in genbank format, with your new feature changed.
print("\nThis bit is the SeqRecord, printed out in genbank format, with a feature changed.\n")
print(my_sequence_record.format("gb"))