from gctox.model import run_inference

if __name__ == '__main__':
    smiles = ['C', 'CC', 'CCC']
    print(run_inference(smiles))

