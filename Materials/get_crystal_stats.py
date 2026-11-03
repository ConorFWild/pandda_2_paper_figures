import joblib
import subprocess
import fire
import pathlib

import gemmi
import pandas as pd
from rich import print as rprint


def get_i_over_sigma(mtz_path):
    mtz = gemmi.read_mtz_file(str(mtz_path))
    binner = gemmi.Binner()
    binner.setup(60, gemmi.Binner.Method.Dstar3, mtz)
    bins = binner.get_bins(mtz)
    fo = mtz.column_with_label('FP').array
    sigma = mtz.column_with_label('SIGFP').array

    means = gemmi.binmean(bins, fo/sigma)
    for idx, mean in enumerate(means):
        if mean < 2:
            dmin = binner.dmin_of_bin(idx-1)

    return dmin

def get_num_asu_res(path):
    pdb = gemmi.read_structure(str(path))
    pdb.setup_entities()
    num_res = 0
    for entity in pdb.entities: 
        num_res += entity.full_sequence
        print(entity.full_sequence)
        print(len(entity.full_sequence))
    return num_res

def get_num_waters(path):
    pdb = gemmi.read_structure(str(path))
    num_waters = 0
    for model in pdb:
        for chain in model:
            for res in chain:
                if res.name == 'HOH':
                    num_waters +=1
    return num_waters 

def get_asu_vol(path):
    pdb = gemmi.read_structure(str(path))
    cell = pdb.cell
    return cell.volume

def main(path, out_path):
    records = []
    for system_dir in pathlib.Path(path).glob('*'):
        if not system_dir.is_dir():
            continue
        print(system_dir.name)
        for dataset_dir in (system_dir / 'processed_datasets').glob('*'):
            if not dataset_dir.is_dir():
                rprint(f'\tNOT A DIR: {dataset_dir}')
                continue
            pdb_file = dataset_dir / f'{dataset_dir.name}-pandda-input.pdb'
            mtz_file = dataset_dir / f'{dataset_dir.name}-pandda-input.mtz'
            if not pdb_file.exists():
                rprint(f"\tNO SUCH FILE: {pdb_file}")
                continue
            if not mtz_file.exists():
                rprint(f"\tNO SUCH FILE: {mtz_file}")
                continue

            i_sigma = get_i_over_sigma(mtz_file)
            asu_res = get_num_asu_res(pdb_file)
            num_waters = get_num_waters(pdb_file)
            asu_vol = get_asu_vol(pdb_file)

            record = {
                'System': system_dir.name,
                'Dtag': dataset_dir.name,
                'I/Sigma=2': i_sigma,
                'ASU Residues': asu_res,
                'ASU Volume': num_waters,
                'ASU Waters': asu_vol,
            }
            rprint(record)
            records.append(record)

    print(f'Finished!')


    df = pd.DataFrame(record)
    df.write_csv(out_path)
    ...

if __name__ == "__main__":
    fire.Fire(main)