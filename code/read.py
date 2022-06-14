import json
import pandas as pd
 
with open("elements.json", "r") as f:
    elements = json.load(f)

df = pd.DataFrame.from_dict(elements)
print(df)
