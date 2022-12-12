import pandas as pd
import plotly.express as px

from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
prof = ParsedParameterFile("4/uniform/profiling")
print(prof["profiling"]) # print python dict

df = pd.DataFrame(prof["profiling"]).T
df.index = df.index.set_names(['trigger'])
df = df.reset_index("trigger")
df.to_csv("performance.csv",index=False)

fig = px.treemap(df,names='id', parents = 'parentId',hover_name='description', values='totalTime',branchvalues="total")
fig.show()

