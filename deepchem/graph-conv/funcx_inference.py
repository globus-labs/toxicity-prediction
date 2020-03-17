"""Run toxicity inference with FuncX"""
from funcx.sdk.client import FuncXClient
import json

# Get FuncX ready
fxc = FuncXClient()
theta_ep = 'd3a23590-3282-429a-8bce-e0ca0f4177f3'
with open('func_uuid.json') as fp:
    func_id = json.load(fp)

# Run the infernece
smiles = ['C', 'CC', 'CCC']
task_id = fxc.run(smiles, endpoint_id=theta_ep, function_id=func_id)
print(task_id)
