
# GISAID Datatframes

GISAID dataframes for use in the CoronaTrend app at [coronatrend.live](coronatrend.live)

Dataframes are generated using _dataframegen.py_.

HOW TO USE:

`python dataframegen.py --aln Aligned.fasta --feather Output.feather`

This outputs the dataframe in .feather format, which can be read by the pandas library in Python.



NOTE 1: 
Sequences where a position has an ambiguous base (e.g. N) are treated as no mutations.

NOTE 2:
Weeks with a total sequence count less than 10 are removed to maintain accuracy of mutation prevalence.


