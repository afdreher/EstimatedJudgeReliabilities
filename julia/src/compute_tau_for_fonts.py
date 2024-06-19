import argparse
import pandas as pd


parser = argparse.ArgumentParser(
    prog='ComputeTaus',
    description='Compute tau between fonts'
)

parser.add_argument('first') 
parser.add_argument('second')

args = parser.parse_args()

first_data = pd.read_csv(args.first, index_col='name')
second_data = pd.read_csv(args.second, index_col='name')

results = {}
for attrib in list(first_data.columns):
    corr_val = pd.DataFrame([first_data[attrib].rename('first'), second_data[attrib].rename('second')]).transpose().corr()['first']['second']
    results[attrib] = {'tau' : corr_val}


dfResults = pd.DataFrame.from_dict(results).transpose()
dfResults.sort_index(inplace=True)

print(dfResults)