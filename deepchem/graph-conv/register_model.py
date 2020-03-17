from funcx.sdk.client import FuncXClient
import json

def wrapped_func(*args):
    from gctox import run_inference
    return run_inference(*args)

fxc = FuncXClient()
func_uuid = fxc.register_function(wrapped_func, description="Infer toxicity based on Tox21 with Deepchem's Graph Convolution")

print(f'Registered function as {func_uuid}')
with open('func_uuid.json', 'w') as fp:
    json.dump(func_uuid, fp)

