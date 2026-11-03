import joblib
import subprocess
import fire
import pathlib

import gemmi
import pandas as pd
from rich import print as rprint

F_LABELS = ['F']
SIG_LABELS = ['SIGF']
RESIDUE_NAMES = ["ALA",
                 "ARG",
                 "ASN",
                 "ASP",
                 "CYS",
                 "GLN",
                 "GLU",
                 "HIS",
                 "ILE",
                 "LEU",
                 "LYS",
                 "MET",
                 "PHE",
                 "PRO",
                 "SER",
                 "THR",
                 "TRP",
                 "TYR",
                 "VAL",
                 "GLY",
                 ]

def get_i_over_sigma(mtz_path):
    mtz = gemmi.read_mtz_file(str(mtz_path))
    binner = gemmi.Binner()
    binner.setup(60, gemmi.Binner.Method.Dstar3, mtz)
    bins = binner.get_bins(mtz)
    fo = None
    for f_label in F_LABELS:
        fo_raw = mtz.column_with_label(f_label)
        if fo_raw is not None:
            fo = fo_raw.array
    if fo is None:
        raise Exception()

    sigma = None
    for sig_label in SIG_LABELS:
        sig_raw = mtz.column_with_label(sig_label)
        if sig_raw is not None:
            sigma = sig_raw.array
    if sigma is None:
        raise Exception()

    means = gemmi.binmean(bins, fo/sigma)
    dmin = None
    for idx, mean in enumerate(means):
        print(f'RES: {round(binner.dmin_of_bin(idx) ,2 )} I/SIGMA {round(mean, 2)}')
        if mean < 2:
            dmin = binner.dmin_of_bin(idx-1)
            break

    if dmin is None:
        dmin = binner.dmin_of_bin(idx)

    return dmin

def get_num_asu_res(path):
    pdb = gemmi.read_structure(str(path))
    pdb.setup_entities()
    num_res = 0
    for model in pdb:
        for chain in model:
            for res in chain:
                if res.name in RESIDUE_NAMES:
                    num_res +=1
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
        print(f'##### {system_dir.name} #####')
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
                'I/Sigma=2': round(i_sigma, 2),
                'ASU Residues': asu_res,
                'ASU Volume': round(asu_vol, 2),
                'ASU Waters': num_waters,
            }
            rprint(record)
            records.append(record)

    print(f'Finished!')


    df = pd.DataFrame(record)
    df.write_csv(out_path)
    ...

if __name__ == "__main__":
    fire.Fire(main)